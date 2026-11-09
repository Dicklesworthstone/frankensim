//! G1/G3: distributed acceleration, fixed support and finite end heat storage.
use super::*;
use crate::body_force::{BodyForce, ForceControl};

const FORCE_CONTROL: ForceControl = ForceControl { maximum_mach: 80.0 };

#[test]
fn g1_heated_reacting_piston_separates_body_work_support_and_finite_storage() {
    let mut flow = closed(&mechanism(8.0, 30.0, 100.0), 4, chemistry());
    let mut operator = molecular(flow.flow());
    let mut thermal = PistonThermalEnds::new(
        flow.flow(),
        [ADIABATIC, finite(3.0, 0.4, 1.0)],
        [None, None],
        || false,
    )
    .unwrap();
    let mut piston = Piston::new(0.2, 0.1).unwrap();
    let before = flow.flow().totals();
    let acceleration = -2.0;
    let report = flow
        .step_with_piston_dissipation_heat_and_body_force(
            &mut operator,
            &mut thermal,
            &mut piston,
            Time::new(0.005),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            WallControl::default(),
            load(),
            piston_controls(),
            BodyForce::UniformAcceleration(acceleration),
            FORCE_CONTROL,
            || false,
        )
        .unwrap();
    let h = report.gas.reactive.elapsed.value();
    let after = flow.flow().totals();
    let body = report.gas.reactive.transport.source_change;
    close(body.impulse, before.mass() * acceleration * h, 3e-12);
    let gas_change = physical_change(flow.flow().gas(), before, after);
    let kinetic = 0.5 * piston.mass_kg() * (piston.velocity_m_s().powi(2) - 0.1_f64.powi(2));
    let thermal_change = 0.4 * (wall_temperature(thermal.contacts()[1]) - 3.0);
    let displacement = flow.flow().mesh().positions()[4] - 1.0;
    close(
        gas_change + kinetic + thermal_change + AREA * displacement - body.work,
        0.0,
        3e-8,
    );
    close(
        after.momentum - before.momentum
            + piston.mass_kg() * (piston.velocity_m_s() - 0.1)
            + report.fixed_wall_impulse_n_s
            + AREA * h
            - body.impulse,
        0.0,
        4e-11,
    );
    close(report.gas.gas_energy_change_j, gas_change, 3e-11);
    close(
        report.total_energy_residual_j,
        0.0,
        piston_controls().energy_tolerance_j,
    );
    assert_eq!(report.gas.heat.external_heat_to_gas_j, 0.0);
    assert_eq!(flow.flow().mesh().positions()[0], 0.0);
    assert_eq!(operator.mesh(), flow.flow().mesh());
    assert_eq!(thermal.mesh(), flow.flow().mesh());
    assert!(report.gas.reactive.chemistry.species_change_kg[1] > 0.0);
}

#[test]
fn g3_zero_source_replays_thermal_piston_and_adiabatic_entry_replays_molecular_force() {
    for forced in [false, true] {
        let mechanism = mechanism(8.0, 30.0, 0.0);
        let mut a = closed(&mechanism, 4, chemistry());
        let mut b = closed(&mechanism, 4, chemistry());
        let mut da = molecular(a.flow());
        let mut db = molecular(b.flow());
        let mut pa = Piston::new(0.2, 0.0).unwrap();
        let mut pb = pa;
        let contacts = if forced {
            [ADIABATIC; 2]
        } else {
            [ADIABATIC, finite(3.0, 0.4, 1.0)]
        };
        let mut ta = PistonThermalEnds::new(a.flow(), contacts, [None, None], || false).unwrap();
        let mut tb = PistonThermalEnds::new(b.flow(), contacts, [None, None], || false).unwrap();
        let force = BodyForce::UniformAcceleration(if forced { 3.0 } else { 0.0 });
        let report = b
            .step_with_piston_dissipation_heat_and_body_force(
                &mut db,
                &mut tb,
                &mut pb,
                Time::new(0.005),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                WallControl::default(),
                load(),
                piston_controls(),
                force,
                FORCE_CONTROL,
                || false,
            )
            .unwrap();
        if forced {
            let direct = a
                .step_with_piston_dissipation_and_body_force(
                    &mut da,
                    &mut pa,
                    Time::new(0.005),
                    ReactiveStepControl::default(),
                    DissipationControl::default(),
                    load(),
                    piston_controls(),
                    force,
                    FORCE_CONTROL,
                    || false,
                )
                .unwrap();
            assert_eq!(report.gas.reactive, direct.gas);
            assert_eq!(report.gas_impulse_n_s, direct.gas_impulse_n_s);
            assert_eq!(
                report.total_energy_residual_j,
                direct.total_energy_residual_j
            );
        } else {
            let direct = step(&mut a, &mut da, &mut ta, &mut pa, || false).unwrap();
            assert_eq!(report, direct);
            assert_eq!(tb.contacts(), ta.contacts());
        }
        assert_eq!(a.accepted.states, b.accepted.states);
        assert_eq!(a.flow().mesh(), b.flow().mesh());
        assert_eq!(pa, pb);
    }
}
