//! G0 bounded work/rollback and G1 automatic coupled endpoint solutions.
use super::*;
use crate::enthalpy::evaporation::implicit::solve::{
    EvaporationImplicitControls, EvaporationImplicitError,
};
use std::cell::Cell;

fn controls() -> EvaporationImplicitControls {
    EvaporationImplicitControls::new(range(0.0))
}

#[test]
fn automatic_stiff_two_cell_solve_matches_the_independent_linear_solution() {
    let dt = 1000.0;
    let mut a = cell(1, 0.1, 300.0, 0.0, 10.0);
    a.boundary.ambient_k = 310.0;
    a.boundary.htc_w_m2_k = 50.0;
    a.boundary.reservoir = Some(Reservoir {
        area_m2: 0.002,
        conductance_w_k: 0.3,
        temperature_k: 320.0,
    });
    let mut b = cell(2, 0.2, 310.0, 0.0, -2.0);
    b.boundary.ambient_k = 295.0;
    b.boundary.htc_w_m2_k = 20.0;
    let old = simulation(&CALORIC, vec![a, b], Some(2.0), vec![]);
    let dm = transfer(&CONSTANT, 300.0, film(0.0), 0.006, dt);
    let diagonal = [400.0 + dt * 2.6, 800.0 + dt * 2.16];
    let rhs = [
        dt * (10.0 + 0.3 * 10.0 + 0.3 * 20.0 + 2.0 * 10.0) - dm * LATENT,
        dt * (-2.0 + 0.16 * -15.0 + 2.0 * -10.0),
    ];
    let determinant = diagonal[0] * diagonal[1] - (dt * 2.0).powi(2);
    let expected = [
        300.0 + (rhs[0] * diagonal[1] + dt * 2.0 * rhs[1]) / determinant,
        310.0 + (diagonal[0] * rhs[1] + dt * 2.0 * rhs[0]) / determinant,
    ];
    assert!(dt > 8.0 * 0.8 * 400.0 / 2.6);
    let calls = Cell::new(0usize);
    let (next, payload, report) = old
        .solve_implicit_evaporating_step(
            CellId(1),
            &CONSTANT,
            film(0.0),
            dt,
            controls(),
            |thermal, liquid| {
                calls.set(calls.get() + 1);
                assert_eq!(thermal.cells()[&CellId(1)].mass_kg, 0.1);
                assert!(liquid.mass_kg < 0.1);
                Ok(coefficients(calls.get()))
            },
            || Ok(()),
        )
        .unwrap();
    close(next.temperature_k(CellId(1)).unwrap(), expected[0], 1e-9);
    close(next.temperature_k(CellId(2)).unwrap(), expected[1], 1e-9);
    close(report.physical.surface.transferred_mass_kg, dm, 1e-15);
    let external = dt
        * (8.0
            + 0.3 * (310.0 - expected[0])
            + 0.3 * (320.0 - expected[0])
            + 0.16 * (295.0 - expected[1]));
    close(report.physical.external_heat_j, external, 1e-6);
    close(
        next.energy_j() - old.energy_j() + report.physical.surface.vapor_energy_j,
        external,
        1e-6,
    );
    assert_eq!(payload, calls.get());
    assert_eq!(report.residual_evaluations, calls.get());
    assert!(report.newton_iterations > 0);
    assert!(report.residual_evaluations <= controls().maximum_evaluations);
    assert_eq!(next.links(), old.links());
    assert_eq!(next.time_s(), dt);
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn nonlinear_mass_dependent_heat_solves_inside_the_liquid_chart_and_preserves_datum() {
    let dt = 1000.0;
    let target = 315.0;
    let dm = transfer(&VARYING, target, film(0.0), 0.01, dt);
    let q = 0.1 * CP * 15.0 + dm * (LATENT - 2000.0 * 15.0);
    assert!(300.0 + q / (0.1 * CP) > 330.0);
    let authored = q / dt + 60.0 * dm;
    let mut receipts = Vec::new();
    for datum in [0.0, 1e9] {
        let caloric = Caloric {
            datum,
            fixed_mass: None,
        };
        let old = simulation(
            &caloric,
            vec![cell(1, 0.1, 300.0, datum, authored)],
            None,
            vec![],
        );
        let mut config = controls();
        config.pure_liquid_enthalpy_range_j_kg = range(datum);
        let (next, (temperature, mass), report) = old
            .solve_implicit_evaporating_step(
                CellId(1),
                &VARYING,
                film(0.0),
                dt,
                config,
                |thermal, liquid| {
                    let mut result =
                        coefficients((thermal.temperature_k(CellId(1))?, liquid.mass_kg));
                    result
                        .transport
                        .heat
                        .bulk
                        .sources_w
                        .insert(CellId(1), 60.0 * (liquid.mass_kg - 0.1));
                    Ok(result)
                },
                || Ok(()),
            )
            .unwrap();
        close(temperature, target, 1e-9);
        close(next.temperature_k(CellId(1)).unwrap(), target, 1e-9);
        close(mass, 0.1 - dm, 1e-12);
        close(report.physical.external_heat_j, q, 1e-6);
        close(
            report.physical.energy_residual_j,
            0.0,
            config.total_energy_tolerance_j,
        );
        assert!(report.newton_iterations >= 2);
        assert_eq!(next.cells()[&CellId(1)].boundary.source_w, authored);
        receipts.push(report.physical);
    }
    close(
        receipts[0].surface.liquid.mass_kg,
        receipts[1].surface.liquid.mass_kg,
        1e-12,
    );
    close(
        receipts[1].surface.vapor_energy_j - receipts[0].surface.vapor_energy_j,
        receipts[0].surface.transferred_mass_kg * 1e9,
        1e-4,
    );
}

#[test]
fn automatic_condensation_uses_the_endpoint_latent_heat_and_increases_inventory() {
    let dt = 1000.0;
    let target = 306.0;
    let latent = LATENT - 2000.0 * 6.0;
    let dm = -0.1 * CP * 6.0 / latent;
    // Invert the independent Stefan-film equation so the no-heater endpoint
    // closes exactly at the chosen temperature, without supplying its h iterate.
    let scale = 0.01 * 0.018 * 100_000.0 / (MOLAR_GAS_CONSTANT * 300.0);
    let pressure = 100_000.0 - (100_000.0 - 4600.0) * (dm / (0.01 * dt * scale)).exp();
    let old = simulation(&CALORIC, vec![cell(1, 0.1, 300.0, 0.0, 0.0)], None, vec![]);
    let (next, (), report) = old
        .solve_implicit_evaporating_step(
            CellId(1),
            &VARYING,
            film(pressure),
            dt,
            controls(),
            empty,
            || Ok(()),
        )
        .unwrap();
    close(next.temperature_k(CellId(1)).unwrap(), target, 1e-9);
    close(report.physical.surface.transferred_mass_kg, dm, 1e-12);
    assert!(next.cells()[&CellId(1)].mass_kg > 0.1);
    assert_eq!(report.physical.external_heat_j, 0.0);
    close(
        next.energy_j() + report.physical.surface.vapor_energy_j,
        0.0,
        1e-6,
    );
}

#[test]
fn work_limits_original_source_refusals_and_physical_gates_survive_newton() {
    let old = simulation(&CALORIC, vec![cell(1, 0.1, 300.0, 0.0, 15.0)], None, vec![]);
    let mut limited = controls();
    limited.maximum_evaluations = 1;
    let calls = Cell::new(0usize);
    let result = old.solve_implicit_evaporating_step(
        CellId(1),
        &VARYING,
        film(0.0),
        1000.0,
        limited,
        |_, _| {
            calls.set(calls.get() + 1);
            Ok(coefficients(()))
        },
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(EvaporationImplicitError::Physical(
            ThermalEvaporationError::Thermal(DepositionError::Budget(
                "implicit evaporation residual evaluations"
            ))
        ))
    ));
    assert_eq!(calls.get(), 1);
    limited = controls();
    limited.maximum_iterations = 1;
    assert!(matches!(
        old.solve_implicit_evaporating_step(
            CellId(1),
            &VARYING,
            film(0.0),
            1000.0,
            limited,
            empty,
            || Ok(()),
        ),
        Err(EvaporationImplicitError::NonConverged(_))
    ));
    let mut loose = controls();
    loose.newton.absolute_tolerance = 1e12;
    assert!(matches!(
        old.solve_implicit_evaporating_step(
            CellId(1),
            &VARYING,
            film(0.0),
            1000.0,
            loose,
            empty,
            || Ok(()),
        ),
        Err(EvaporationImplicitError::Physical(
            ThermalEvaporationError::Thermal(DepositionError::EnergyBalance(_))
        ))
    ));
    let calls = Cell::new(0usize);
    let result = old.solve_implicit_evaporating_step(
        CellId(1),
        &VARYING,
        film(0.0),
        1000.0,
        controls(),
        |thermal, _| {
            calls.set(calls.get() + 1);
            if thermal.temperature_k(CellId(1))? > 300.0 {
                Err(ThermalEvaporationError::Surface(EvaporationError::Invalid(
                    "original source stencil refused",
                )))
            } else {
                Ok(coefficients(()))
            }
        },
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(EvaporationImplicitError::Physical(
            ThermalEvaporationError::Surface(EvaporationError::Invalid(
                "original source stencil refused"
            ))
        ))
    ));
    assert_eq!(calls.get(), 2);
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.cells()[&CellId(1)].mass_kg, 0.1);
    assert_eq!(old.cells()[&CellId(1)].specific_h_j_kg, 0.0);
}

#[test]
fn final_physical_recheck_is_fresh_funded_and_cancellable() {
    let source = transfer(&CONSTANT, 300.0, film(0.0), 0.01, 1.0) * LATENT;
    let old = simulation(
        &CALORIC,
        vec![cell(1, 0.1, 300.0, 0.0, source)],
        None,
        vec![],
    );
    let calls = Cell::new(0usize);
    let result = old.solve_implicit_evaporating_step(
        CellId(1),
        &CONSTANT,
        film(0.0),
        1.0,
        controls(),
        |_, _| {
            calls.set(calls.get() + 1);
            if calls.get() == 2 {
                Err(ThermalEvaporationError::Surface(EvaporationError::Invalid(
                    "final source recheck refused",
                )))
            } else {
                Ok(coefficients(()))
            }
        },
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(EvaporationImplicitError::Physical(
            ThermalEvaporationError::Surface(EvaporationError::Invalid(
                "final source recheck refused"
            ))
        ))
    ));
    assert_eq!(calls.get(), 2);
    let mut limited = controls();
    limited.maximum_evaluations = 1;
    assert!(matches!(
        old.solve_implicit_evaporating_step(
            CellId(1),
            &CONSTANT,
            film(0.0),
            1.0,
            limited,
            empty,
            || Ok(()),
        ),
        Err(EvaporationImplicitError::Physical(
            ThermalEvaporationError::Thermal(DepositionError::Budget(
                "implicit evaporation residual evaluations"
            ))
        ))
    ));
    let mut exact_mass = controls();
    exact_mass.mass_tolerance_kg = 0.0;
    assert!(matches!(
        old.solve_implicit_evaporating_step(
            CellId(1),
            &CONSTANT,
            film(0.0),
            1.0,
            exact_mass,
            empty,
            || Ok(()),
        ),
        Err(EvaporationImplicitError::Physical(
            ThermalEvaporationError::MassBalance { .. }
        ))
    ));

    let old = simulation(&CALORIC, vec![cell(1, 0.1, 300.0, 0.0, 15.0)], None, vec![]);
    let mut polls = 0;
    old.solve_implicit_evaporating_step(
        CellId(1),
        &VARYING,
        film(0.0),
        1000.0,
        controls(),
        empty,
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    for cancel_at in [1, polls / 2, polls] {
        let mut count = 0;
        let result = old.solve_implicit_evaporating_step(
            CellId(1),
            &VARYING,
            film(0.0),
            1000.0,
            controls(),
            empty,
            || {
                count += 1;
                if count == cancel_at {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        );
        assert!(matches!(
            result,
            Err(EvaporationImplicitError::Physical(
                ThermalEvaporationError::Thermal(DepositionError::Cancelled)
            ))
        ));
    }
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.cells()[&CellId(1)].mass_kg, 0.1);
}
