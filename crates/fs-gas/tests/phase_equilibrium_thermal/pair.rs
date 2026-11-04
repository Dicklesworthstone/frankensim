//! Two finite source inventories must exchange heat without replacing either law.

use super::{
    IsobaricMixtureEquilibriumParcel, MIXTURE_MASS, Pressure, Temperature,
    ThermochemicalPengRobinsonModel, Time, WallExchangeError, close, equilibrium,
    isobaric_controls, isobaric_parcel, model, parcel, source,
};
use fs_gas::material::phase_thermal::{
    IsobaricMixturePairThermalControl, MixedMixturePairThermalControl, MixturePairThermalControl,
    exchange_isobaric_mixture_pair, exchange_isochoric_isobaric_mixture_pair,
    exchange_isochoric_mixture_pair,
};

fn controls(bracket: [f64; 2]) -> IsobaricMixturePairThermalControl {
    let phase = isobaric_controls(bracket).equilibrium;
    IsobaricMixturePairThermalControl::new(bracket, phase, phase)
}

fn at_pressure(
    model: ThermochemicalPengRobinsonModel,
    mass: f64,
    temperature: f64,
    pressure: f64,
) -> IsobaricMixtureEquilibriumParcel {
    IsobaricMixtureEquilibriumParcel::from_pressure_temperature(
        model,
        mass,
        Temperature::new(temperature),
        Pressure::new(pressure),
        controls([-1.0, 1.0]).left_equilibrium.phase_options(),
        || false,
    )
    .unwrap()
}

#[test]
fn g1_finite_phase_pair_closes_heat_and_both_pressure_work_accounts() {
    let mass = 2.0 * MIXTURE_MASS;
    let mut cold = isobaric_parcel(source(), mass, 180.0);
    let mut hot = at_pressure(source(), mass, 280.0, 2e6);
    let before = [cold.clone(), hot.clone()];
    // Independent binary-mixture values at 3 MPa select a cold endpoint of
    // 220 K, across the original 205.7464 K bubble point.
    let heat = mass * (-16_172.250_165_027_297 + 203_472.419_617_408_38);
    let cold_stored = mass * (-39_099.945_545_128_86 + 209_511.548_518_424_48);
    // The second original parcel supplies an independent material-only
    // reference at its own pressure; it does not use the pair heat solver.
    let mut hot_reference = hot.clone();
    let reference = hot_reference
        .add_heat(
            -heat,
            controls([0.0, 1.01 * heat]).right_equilibrium,
            || false,
        )
        .unwrap();
    assert!(hot_reference.state().temperature_k() > 220.0);
    let duration = Time::new(heat / (hot_reference.state().temperature_k() - 220.0));
    let report = exchange_isobaric_mixture_pair(
        &mut cold,
        &mut hot,
        1.0,
        duration,
        controls([0.0, 1.01 * heat]),
        || false,
    )
    .unwrap();
    close(cold.state().temperature_k(), 220.0, 2e-6);
    assert!(cold.state().phase().two_phase().is_some());
    close(
        hot.state().temperature_k(),
        hot_reference.state().temperature_k(),
        2e-6,
    );
    close(report.heat_to_left_j, heat, 2e-5);
    close(report.heat_to_right_j, -heat, 2e-5);
    close(report.internal_energy_change_j[0], cold_stored, 2e-5);
    close(
        report.internal_energy_change_j[1],
        reference.stored_energy_change_j,
        2e-5,
    );
    for (index, current) in [&cold, &hot].into_iter().enumerate() {
        let actual_work = current.pressure_pa() * (current.volume_m3() - before[index].volume_m3());
        assert_eq!(report.boundary_work_by_parcels_j[index], actual_work);
        assert_eq!(current.pressure_pa(), before[index].pressure_pa());
        assert_eq!(current.mass_kg(), before[index].mass_kg());
        assert_eq!(current.model(), before[index].model());
    }
    assert!(report.boundary_work_by_parcels_j[0] > 0.0);
    assert!(report.boundary_work_by_parcels_j[1] < 0.0);
    close(
        report.left_equilibrium.enthalpy_change_j + report.right_equilibrium.enthalpy_change_j,
        0.0,
        2e-5,
    );
    close(
        report.internal_energy_change_j.iter().sum::<f64>()
            + report.boundary_work_by_parcels_j.iter().sum::<f64>(),
        0.0,
        2e-5,
    );
    close(report.energy_residual_j, 0.0, 2e-5);
    close(report.thermal_residual_j, 0.0, 2e-5);
    assert!(report.entropy_production_j_k > 0.0);
    assert_eq!(
        report.entropy_production_j_k,
        report.entropy_change_j_k.iter().sum::<f64>()
    );
    assert_eq!(report.pressure_evaluations, 0);
    assert!(report.pt_evaluations <= controls([0.0, heat]).maximum_pt_evaluations);
}

#[test]
fn g1_mixed_and_rigid_pairs_retain_each_original_storage_constraint() {
    let mass = 2.0 * MIXTURE_MASS;
    let mut rigid = parcel(source(), mass);
    let mut expanding = isobaric_parcel(source(), mass, 280.0);
    let initial_rigid_volume = rigid.volume_m3();
    let mut expected_rigid = rigid.clone();
    let mut expected_expanding = expanding.clone();
    let heat = 10.0;
    expected_rigid
        .add_heat(heat, equilibrium(), || false)
        .unwrap();
    expected_expanding
        .add_heat(-heat, controls([0.0, 11.0]).right_equilibrium, || false)
        .unwrap();
    let duration = Time::new(
        heat / (expected_expanding.state().temperature_k()
            - expected_rigid.state().temperature_k()),
    );
    let control = MixedMixturePairThermalControl::new(
        [0.0, 11.0],
        equilibrium(),
        controls([0.0, 11.0]).right_equilibrium,
    );
    let report = exchange_isochoric_isobaric_mixture_pair(
        &mut rigid,
        &mut expanding,
        1.0,
        duration,
        control,
        || false,
    )
    .unwrap();
    close(
        rigid.state().temperature_k(),
        expected_rigid.state().temperature_k(),
        2e-6,
    );
    close(
        expanding.state().temperature_k(),
        expected_expanding.state().temperature_k(),
        2e-6,
    );
    assert_eq!(rigid.volume_m3(), initial_rigid_volume);
    assert_eq!(report.boundary_work_by_parcels_j[0], 0.0);
    assert!(report.boundary_work_by_parcels_j[1] < 0.0);
    close(report.left_equilibrium.stored_energy_change_j, heat, 2e-5);
    close(report.right_equilibrium.enthalpy_change_j, -heat, 2e-5);
    close(report.energy_residual_j, 0.0, 2e-5);
    assert!(report.entropy_production_j_k > 0.0);

    // A second rigid parcel uses the same kernel with both work accounts zero.
    let mut left = parcel(source(), mass);
    let mut right = parcel(source(), mass);
    right.add_heat(50.0, equilibrium(), || false).unwrap();
    let volumes = [left.volume_m3(), right.volume_m3()];
    let rigid_report = exchange_isochoric_mixture_pair(
        &mut left,
        &mut right,
        1.0,
        Time::new(0.01),
        MixturePairThermalControl::new([0.0, 1.0], equilibrium(), equilibrium()),
        || false,
    )
    .unwrap();
    assert_eq!([left.volume_m3(), right.volume_m3()], volumes);
    assert_eq!(rigid_report.boundary_work_by_parcels_j, [0.0; 2]);
    assert!(rigid_report.heat_to_left_j > 0.0);
    close(
        rigid_report.internal_energy_change_j.iter().sum(),
        0.0,
        2e-5,
    );
}

#[test]
fn g3_pair_heat_and_entropy_ignore_independent_species_reference_datums() {
    let mass = 2.0 * MIXTURE_MASS;
    let mut left = isobaric_parcel(source(), mass, 220.0);
    let mut right = isobaric_parcel(source(), mass, 280.0);
    let mut shifted_left = isobaric_parcel(
        model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]),
        mass,
        220.0,
    );
    let mut shifted_right = isobaric_parcel(
        model([-4e21, 9e21], [-8e17, 3e18], [100.0, 450.0]),
        mass,
        280.0,
    );
    let control = controls([0.0, 0.1]);
    let a = exchange_isobaric_mixture_pair(
        &mut left,
        &mut right,
        1.0,
        Time::new(1e-3),
        control,
        || false,
    )
    .unwrap();
    let b = exchange_isobaric_mixture_pair(
        &mut shifted_left,
        &mut shifted_right,
        1.0,
        Time::new(1e-3),
        control,
        || false,
    )
    .unwrap();
    assert_eq!(a, b);
    assert_eq!(
        left.state().temperature_k(),
        shifted_left.state().temperature_k()
    );
    assert_eq!(
        right.state().temperature_k(),
        shifted_right.state().temperature_k()
    );
    assert_eq!(left.volume_m3(), shifted_left.volume_m3());
    assert_eq!(right.volume_m3(), shifted_right.volume_m3());
    assert!(a.entropy_production_j_k > 0.0);
}

#[test]
fn g4_pair_rest_refusals_and_final_cancellation_preserve_both_inventories() {
    let mass = 2.0 * MIXTURE_MASS;
    let original = [
        isobaric_parcel(source(), mass, 220.0),
        isobaric_parcel(source(), mass, 280.0),
    ];
    let [mut left, mut right] = original.clone();
    let normal = controls([0.0, 100.0]);
    let mut limited = normal;
    limited.maximum_pt_evaluations = 1;
    let rest =
        exchange_isobaric_mixture_pair(&mut left, &mut right, 1.0, Time::new(0.0), limited, || {
            false
        })
        .unwrap();
    assert_eq!(rest.pt_evaluations, 0);
    assert_eq!(rest.heat_to_left_j, 0.0);
    assert_eq!([left.clone(), right.clone()], original);
    assert!(matches!(
        exchange_isobaric_mixture_pair(&mut left, &mut right, 1.0, Time::new(1.0), limited, || {
            false
        }),
        Err(WallExchangeError::MixtureThermalWork { .. })
    ));
    assert_eq!([left.clone(), right.clone()], original);
    let mut narrow = normal;
    narrow.heat_bracket_j = [0.0, 0.001];
    assert!(
        exchange_isobaric_mixture_pair(&mut left, &mut right, 1.0, Time::new(1.0), narrow, || {
            false
        })
        .is_err()
    );
    assert_eq!([left.clone(), right.clone()], original);
    let mut polls = 0;
    let expected =
        exchange_isobaric_mixture_pair(&mut left, &mut right, 1.0, Time::new(1.0), normal, || {
            polls += 1;
            false
        })
        .unwrap();
    let [mut retry_left, mut retry_right] = original.clone();
    let mut seen = 0;
    assert_eq!(
        exchange_isobaric_mixture_pair(
            &mut retry_left,
            &mut retry_right,
            1.0,
            Time::new(1.0),
            normal,
            || {
                seen += 1;
                seen == polls
            },
        ),
        Err(WallExchangeError::Cancelled)
    );
    assert_eq!([retry_left.clone(), retry_right.clone()], original);
    let retried = exchange_isobaric_mixture_pair(
        &mut retry_left,
        &mut retry_right,
        1.0,
        Time::new(1.0),
        normal,
        || false,
    )
    .unwrap();
    assert_eq!(retried, expected);
    assert_eq!([retry_left, retry_right], [left, right]);
}
