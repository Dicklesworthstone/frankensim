//! Synthetic source mixtures: independent PR/caloric and shared-wall references.

use super::*;
use crate::cavity_phs::piston::tests::real_gas::{real_controls as single_controls, real_parcel};
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
    energy::parcel::IsochoricRealGasParcel,
};
use fs_ports::PortedModule;
use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseV1;
use fs_thermochem::{Composition, CompositionBasis};

fn parcel(x_a: f64, temperature_k: f64, pressure_pa: f64) -> IsochoricRealGasParcel {
    let original = real_parcel(temperature_k, pressure_pa);
    let source = original.model().source();
    let mechanics = original.model().mechanics();
    let bindings = source
        .thermodynamics()
        .components()
        .iter()
        .zip(mechanics.components())
        .map(|(source, critical)| PengRobinsonSpeciesBinding {
            species: source.species().clone(),
            critical_temperature: critical.critical_temperature,
            critical_pressure: critical.critical_pressure,
            acentric_factor: critical.acentric_factor,
        })
        .collect();
    let source = source
        .with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![x_a, 1.0 - x_a]).unwrap(),
        )
        .unwrap();
    let model = ThermochemicalPengRobinsonModel::new(
        source,
        bindings,
        mechanics.binary_interactions().to_vec(),
    )
    .unwrap();
    IsochoricRealGasParcel::new(
        model,
        original.point().clone(),
        original.volume_m3(),
        PengRobinsonPhaseV1::Vapor,
    )
    .unwrap()
}

pub(in crate::cavity_phs::piston) fn real_pair() -> DoubleActingGasPiston<IsochoricRealGasParcel> {
    DoubleActingGasPiston::new(
        GasPiston::new(real_parcel(300.0, 5e6), 2.0, 0.01, 0.4).unwrap(),
        parcel(0.35, 350.0, 3.8e6),
        0.008,
    )
    .unwrap()
}

pub(in crate::cavity_phs::piston) fn real_controls() -> DoublePistonSolve {
    DoublePistonSolve {
        piston: single_controls(),
        opposing_gas: single_controls().gas,
        volume_tolerance_m3: 1e-18,
    }
}

fn input(heat_j: [f64; 2]) -> DoublePistonInput {
    DoublePistonInput {
        duration_s: 0.00002,
        heat_j,
        external_force_n: -100.0,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}"
    );
}

fn energy_change(old: &IsochoricRealGasParcel, new: &IsochoricRealGasParcel) -> f64 {
    old.mass_kg()
        * old
            .model()
            .specific_internal_energy_change_between_densities(
                old.density_kg_m3(),
                old.state().temperature_k(),
                new.density_kg_m3(),
                new.state().temperature_k(),
            )
            .unwrap()
}

#[test]
fn dense_chambers_match_independent_joint_caloric_and_impulse_roots() {
    let mut wall = real_pair();
    let initial = wall.clone();
    let drive = input([80.0, -30.0]);
    let report = wall.step(drive, real_controls(), || false).unwrap();
    // Independent 80-digit Decimal solve: both PR/NASA midpoint energy
    // equations and M*dv=dt*(A1*pm1-A2*pm2+F), not a fixed-pressure load.
    // Each caloric endpoint separately obeys u1-u0=2*(um-u0).
    near(report.displacement_m, 9.973_706_725_084_38e-6, 3e-13);
    near(wall.primary().velocity_m_s(), 0.597_370_672_508_438_1, 3e-9);
    let expected = [
        (
            0.001_000_099_737_067_250_8,
            300.864_029_514_699_9,
            5_017_507.639_051_746,
            301.727_871_360_385_9,
            5_035_006.790_219_382,
            0.500_431_496_827_726_5,
            79.499_568_503_172_27,
        ),
        (
            0.000_999_920_210_346_199_3,
            349.332_640_522_623_5,
            3_792_251.142_459_206_8,
            348.665_213_625_305_5,
            3_784_500.414_912_209_5,
            -0.302_582_405_782_034_5,
            -29.697_417_594_217_965,
        ),
    ];
    let old = [initial.primary().gas(), initial.opposing()];
    let new = [wall.primary().gas(), wall.opposing()];
    for i in 0..2 {
        let (volume, tm, pm, temperature, pressure, work, stored) = expected[i];
        near(new[i].volume_m3(), volume, 3e-15);
        near(report.midpoint[i].temperature_k(), tm, 3e-8);
        near(report.midpoint[i].pressure_pa(), pm, 7e-4);
        near(new[i].state().temperature_k(), temperature, 6e-8);
        near(new[i].state().pressure_pa(), pressure, 1e-3);
        near(report.gas[i].work_by_gas_j, work, 2e-7);
        near(report.gas[i].stored_energy_change_j, stored, 2e-7);
        near(energy_change(old[i], new[i]), stored, 2e-7);
        near(
            report.gas[i].work_by_gas_j,
            report.midpoint[i].pressure_pa() * (new[i].volume_m3() - old[i].volume_m3()),
            2e-10,
        );
        near(
            report.midpoint[i].density_kg_m3(),
            old[i].mass_kg() / f64::midpoint(old[i].volume_m3(), new[i].volume_m3()),
            2e-10,
        );
        assert_eq!(
            new[i].density_kg_m3(),
            new[i].mass_kg() / new[i].volume_m3()
        );
        assert_eq!(new[i].mass_kg(), old[i].mass_kg());
        assert_eq!(new[i].model(), old[i].model());
        assert_eq!(new[i].point().axes().get("specimen-lot"), Some(&42.0));
        assert_ne!(
            report.midpoint[i].density_kg_m3(),
            report.midpoint[i].source().state().density
        );
    }
    assert_ne!(old[0].model(), old[1].model());
    assert_ne!(old[0].density_kg_m3(), old[1].density_kg_m3());
    near(
        new[0].volume_m3() - old[0].volume_m3(),
        initial.primary().area_m2() * report.displacement_m,
        1e-18,
    );
    near(
        new[1].volume_m3() - old[1].volume_m3(),
        -initial.opposing_area_m2() * report.displacement_m,
        1e-18,
    );
    let force = initial.primary().area_m2() * report.midpoint[0].pressure_pa()
        - initial.opposing_area_m2() * report.midpoint[1].pressure_pa()
        + drive.external_force_n;
    near(
        initial.primary().mass_kg() * (report.velocity_m_s[1] - report.velocity_m_s[0]),
        drive.duration_s * force,
        real_controls().piston.momentum_tolerance_n_s,
    );
    let kinetic = 0.5
        * initial.primary().mass_kg()
        * (report.velocity_m_s[1] - report.velocity_m_s[0])
        * (report.velocity_m_s[1] + report.velocity_m_s[0]);
    near(
        energy_change(old[0], new[0]) + energy_change(old[1], new[1]) + kinetic,
        drive.heat_j.iter().sum::<f64>() + drive.external_force_n * report.displacement_m,
        real_controls().piston.energy_tolerance_j,
    );
    assert!(report.total_energy_residual_j.abs() <= real_controls().piston.energy_tolerance_j);
    assert!(report.mechanical_energy_residual_j.abs() <= real_controls().piston.energy_tolerance_j);
    assert!(report.volume_exchange_residual_m3.abs() <= real_controls().volume_tolerance_m3);
    assert!(report.material_trials <= 2 * (real_controls().piston.maximum_iterations + 3));
    wall.audit_last_window().unwrap();
}

#[test]
fn opposing_source_changes_motion_and_its_native_tangent_enters_the_joint_root() {
    let original = real_pair();
    let mut changed = DoubleActingGasPiston::new(
        original.primary().clone(),
        parcel(0.75, 350.0, 3.8e6),
        original.opposing_area_m2(),
    )
    .unwrap();
    let mut baseline = original.clone();
    baseline
        .step(input([0.0, 30.0]), real_controls(), || false)
        .unwrap();
    changed
        .step(input([0.0, 30.0]), real_controls(), || false)
        .unwrap();
    assert_eq!(original.opposing().state().pressure_pa(), 3.8e6);
    assert_ne!(
        baseline.primary().velocity_m_s(),
        changed.primary().velocity_m_s()
    );
    assert_ne!(
        baseline.opposing().state().temperature_k(),
        changed.opposing().state().temperature_k()
    );
    assert_ne!(
        baseline.opposing().state().sound_speed_m_s(),
        changed.opposing().state().sound_speed_m_s()
    );

    let volume = 0.001_000_2;
    let delta = 1e-8;
    let center = original
        .trial(volume, input([80.0, -30.0]), real_controls(), &mut || false)
        .unwrap();
    let plus = original
        .trial(
            volume + delta,
            input([80.0, -30.0]),
            real_controls(),
            &mut || false,
        )
        .unwrap();
    let minus = original
        .trial(
            volume - delta,
            input([80.0, -30.0]),
            real_controls(),
            &mut || false,
        )
        .unwrap();
    let finite_difference = (plus.residual - minus.residual) / (2.0 * delta);
    assert!((center.tangent / finite_difference - 1.0).abs() < 2e-6);
    assert!(center.tangent > center.value.primary.tangent);
}

#[test]
fn dense_pair_reverses_and_continues_with_two_true_density_energy_increments() {
    let initial = real_pair();
    let mut forward = initial.clone();
    forward
        .step(input([80.0, -30.0]), real_controls(), || false)
        .unwrap();
    let mut reverse = DoubleActingGasPiston::new(
        GasPiston::new(
            forward.primary().gas().clone(),
            forward.primary().mass_kg(),
            forward.primary().area_m2(),
            -forward.primary().velocity_m_s(),
        )
        .unwrap(),
        forward.opposing().clone(),
        forward.opposing_area_m2(),
    )
    .unwrap();
    reverse
        .step(input([-80.0, 30.0]), real_controls(), || false)
        .unwrap();
    near(
        reverse.primary().velocity_m_s(),
        -initial.primary().velocity_m_s(),
        1e-8,
    );
    for (old, returned) in [
        (initial.primary().gas(), reverse.primary().gas()),
        (initial.opposing(), reverse.opposing()),
    ] {
        near(returned.volume_m3(), old.volume_m3(), 5e-15);
        near(
            returned.state().temperature_k(),
            old.state().temperature_k(),
            1e-7,
        );
        near(
            returned.state().pressure_pa(),
            old.state().pressure_pa(),
            2e-3,
        );
        assert_eq!(returned.mass_kg(), old.mass_kg());
        assert_eq!(returned.model(), old.model());
    }
    let mut heat = 50.0;
    for _ in 0..15 {
        let drive = input([0.1, -0.05]);
        forward.step(drive, real_controls(), || false).unwrap();
        heat += drive.heat_j.iter().sum::<f64>();
    }
    let displacement = (forward.primary().gas().volume_m3() - initial.primary().gas().volume_m3())
        / initial.primary().area_m2();
    let kinetic = 0.5
        * initial.primary().mass_kg()
        * (forward.primary().velocity_m_s() - initial.primary().velocity_m_s())
        * (forward.primary().velocity_m_s() + initial.primary().velocity_m_s());
    near(
        energy_change(initial.primary().gas(), forward.primary().gas())
            + energy_change(initial.opposing(), forward.opposing())
            + kinetic,
        heat + input([0.0; 2]).external_force_n * displacement,
        16.0 * real_controls().piston.energy_tolerance_j,
    );
    forward.audit_last_window().unwrap();
}

#[test]
fn dense_pair_refusals_and_last_cancellation_preserve_both_materials_and_port_history() {
    let mut original = real_pair();
    original
        .step(input([80.0, -30.0]), real_controls(), || false)
        .unwrap();
    for drive in [
        input([1e9, 0.0]),
        input([0.0, 1e9]),
        DoublePistonInput {
            duration_s: 0.0,
            ..input([0.0; 2])
        },
    ] {
        let mut refused = original.clone();
        assert!(refused.step(drive, real_controls(), || false).is_err());
        assert_eq!(refused, original);
        assert_eq!(
            refused.last_window_transfers(),
            original.last_window_transfers()
        );
    }
    let unsupported = DoublePistonSolve {
        opposing_gas: GasEnergySolve {
            maximum_temperature_k: 340.0,
            ..single_controls().gas
        },
        ..real_controls()
    };
    let mut refused = original.clone();
    assert!(
        refused
            .step(input([1.0, 2.0]), unsupported, || false)
            .is_err()
    );
    assert_eq!(refused, original);
    assert_eq!(
        refused.last_window_transfers(),
        original.last_window_transfers()
    );

    let drive = input([1.0, 2.0]);
    let mut accepted = original.clone();
    let mut polls = 0;
    let expected = accepted
        .step(drive, real_controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    for cancel_at in [1, polls / 2, polls - 1, polls] {
        let mut refused = original.clone();
        let mut count = 0;
        assert!(matches!(
            refused.step(drive, real_controls(), || {
                count += 1;
                count == cancel_at
            }),
            Err(PistonError::Cancelled)
        ));
        assert_eq!(refused, original);
        assert_eq!(
            refused.last_window_transfers(),
            original.last_window_transfers()
        );
        assert_eq!(
            refused.last_window_duration_s(),
            original.last_window_duration_s()
        );
        assert_eq!(
            refused.step(drive, real_controls(), || false).unwrap(),
            expected
        );
        assert_eq!(refused, accepted);
    }
}

#[test]
fn opposing_dense_stiffness_refuses_a_clock_that_every_primary_state_admits() {
    let mut wall = DoubleActingGasPiston::new(
        GasPiston::new(real_parcel(300.0, 1e5), 2.0, 0.01, 0.0).unwrap(),
        parcel(0.35, 350.0, 5e6),
        0.03,
    )
    .unwrap();
    let initial = wall.clone();
    let primary = wall.primary().gas();
    let stiffness = |gas: &IsochoricRealGasParcel, area: f64| {
        gas.state().density_kg_m3() * gas.state().sound_speed_m_s().powi(2) * area.powi(2)
            / gas.volume_m3()
    };
    let primary_rate =
        (stiffness(primary, wall.primary().area_m2()) / wall.primary().mass_kg()).sqrt();
    let coupled_rate = ((stiffness(primary, wall.primary().area_m2())
        + stiffness(wall.opposing(), wall.opposing_area_m2()))
        / wall.primary().mass_kg())
    .sqrt();
    let duration = 0.035 / primary_rate;
    assert!(coupled_rate * duration > fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP);
    let stationary = primary
        .midpoint_heat_volume(primary.volume_m3(), 0.0, single_controls().gas, &mut || {
            false
        })
        .unwrap();
    for state in [
        primary.state(),
        stationary.midpoint(),
        stationary.endpoint().state(),
    ] {
        wall.primary()
            .admit_linear_clock(state, primary.volume_m3(), duration)
            .unwrap();
    }
    let balanced = DoublePistonInput {
        duration_s: duration,
        heat_j: [0.0; 2],
        external_force_n: wall.opposing_area_m2() * wall.opposing().state().pressure_pa()
            - wall.primary().area_m2() * primary.state().pressure_pa(),
    };
    assert!(matches!(
        wall.step(balanced, real_controls(), || false),
        Err(PistonError::TimeResolution(_))
    ));
    assert_eq!(wall, initial);
    assert_eq!(
        wall.last_window_transfers(),
        initial.last_window_transfers()
    );
    let report = wall
        .step(
            DoublePistonInput {
                duration_s: 0.02 / coupled_rate,
                ..balanced
            },
            real_controls(),
            || false,
        )
        .unwrap();
    near(report.displacement_m, 0.0, 1e-15);
    near(wall.primary().velocity_m_s(), 0.0, 1e-10);
}
