//! A free-sided finite body on the original conservative caloric timeline.
//! Independent Hencky equilibrium, quadratic recovery/creep and caloric roots
//! exercise physical heating and reversible thermal storage together.

use super::*;
use fs_material::plastic::finite::thermal::source::recovery::{
    THERMAL_FINITE_RECOVERING_NORTON_J2_LAW_ID, THERMAL_FINITE_RECOVERING_NORTON_LAW_ID,
};

const G_SLOPE: f64 = -40e3;
const K_SLOPE: f64 = -60e3;
const EXPANSION: f64 = 1e-4;
const EXPANSION_SLOPE: f64 = 1e-7;
const PRIMARY_SCALE: f64 = 0.04;
const RECOVERY_RATE: f64 = 0.15;
const RECOVERY_ACTIVATION: f64 = 3000.0;
const YIELD: f64 = 2e6;
const HARDENING: f64 = 40e6;

fn recovering_fixture(combined: bool) -> Fixture {
    let mut source = creep_model(1.0);
    source.law = LawId(
        if combined {
            THERMAL_FINITE_RECOVERING_NORTON_J2_LAW_ID
        } else {
            THERMAL_FINITE_RECOVERING_NORTON_LAW_ID
        }
        .into(),
    );
    source.state_schema_version = 3;
    for (name, value, dims) in [
        ("shear_modulus_slope_pa_k", G_SLOPE, PA_K),
        ("bulk_modulus_slope_pa_k", K_SLOPE, PA_K),
        (
            "logarithmic_expansion_per_k",
            EXPANSION,
            Dims([0, 0, 0, -1, 0, 0]),
        ),
        (
            "logarithmic_expansion_slope_per_k2",
            EXPANSION_SLOPE,
            Dims([0, 0, 0, -2, 0, 0]),
        ),
        (
            "creep_hardening_reference_strain",
            PRIMARY_SCALE,
            Dims::NONE,
        ),
        ("creep_hardening_exponent", 1.0, Dims::NONE),
        ("creep_residual_rate_fraction", 0.0, Dims::NONE),
        (
            "creep_reference_recovery_rate_s",
            RECOVERY_RATE,
            Dims([0, 0, -1, 0, 0, 0]),
        ),
        (
            "creep_recovery_activation_temperature_k",
            RECOVERY_ACTIVATION,
            K,
        ),
    ] {
        source
            .parameters
            .insert(name.into(), LawParameter { value, dims });
    }
    if combined {
        for (name, value, dims) in [
            ("yield_stress_pa", YIELD, PA),
            ("hardening_modulus_pa", HARDENING, PA),
            ("yield_stress_slope_pa_k", 0.0, PA_K),
            ("hardening_modulus_slope_pa_k", 0.0, PA_K),
            ("maximum_accumulated_creep_strain", 3.0, Dims::NONE),
        ] {
            source
                .parameters
                .insert(name.into(), LawParameter { value, dims });
        }
    }
    Fixture::with_model(CP, 19, scalar(DENSITY_PROPERTY, Density::DIMS, RHO), source)
}

#[derive(Clone, Default)]
struct IntervalHeat {
    watts: f64,
    delivered_j: f64,
}
impl DepositionObserver for IntervalHeat {
    type Error = DepositionError;
    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        checkpoint()?;
        Ok([(CELL, self.watts)].into())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        duration_s: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        self.delivered_j += self.watts * duration_s;
        Ok(())
    }
}

#[derive(Clone, Copy, Default)]
struct Reference {
    temperature: f64,
    q: f64,
    a: f64,
    b: f64,
    creep: f64,
    plastic: f64,
    resistance: f64,
    log_thermal_stretch: f64,
    free: f64,
    entropy: f64,
    plastic_heat: f64,
    creep_heat: f64,
    coefficient_change: f64,
    ts_change: f64,
    reversible: f64,
    work: f64,
    defect: f64,
}

fn mechanical_reference(
    old: Reference,
    temperature: f64,
    dt: f64,
    load: f64,
    combined: bool,
) -> Reference {
    let delta_t = temperature - 300.0;
    let g = MU + G_SLOPE * delta_t;
    let k = BULK + K_SLOPE * delta_t;
    let log_theta = EXPANSION * delta_t + 0.5 * EXPANSION_SLOPE * delta_t * delta_t;
    let viscosity = VISCOSITY * (ACTIVATION * (1.0 / temperature - 1.0 / 300.0)).exp();
    let rate = RECOVERY_RATE * (RECOVERY_ACTIVATION * (1.0 / 300.0 - 1.0 / temperature)).exp();
    let divisor = 1.0 + rate * dt;
    let scale = PRIMARY_SCALE * divisor;
    let creep_increment = |q: f64| {
        let shifted = scale + old.resistance;
        let drive = scale * dt * q / viscosity;
        2.0 * drive / (shifted + (shifted * shifted + 4.0 * drive).sqrt())
    };
    let plastic_increment = |q: f64| {
        if combined {
            ((q - YIELD) / HARDENING - old.plastic).max(0.0)
        } else {
            0.0
        }
    };
    let compliance = 1.0 / (9.0 * k) + 1.0 / (3.0 * g);
    let q = if load == 0.0 {
        0.0
    } else {
        let residual = |q: f64| {
            (q / load).ln()
                - log_theta
                - old.creep
                - creep_increment(q)
                - old.plastic
                - plastic_increment(q)
                - compliance * q
        };
        let (mut low, mut high) = (load, 2.0 * load);
        assert!(residual(low) < 0.0 && residual(high) > 0.0);
        for _ in 0..100 {
            let q = 0.5 * (low + high);
            if residual(q) < 0.0 {
                low = q;
            } else {
                high = q;
            }
        }
        0.5 * (low + high)
    };
    let dc = creep_increment(q);
    let dp = plastic_increment(q);
    let creep = old.creep + dc;
    let plastic = old.plastic + dp;
    let a = log_theta + creep + plastic + compliance * q;
    let b = log_theta - 0.5 * (creep + plastic) + q / (9.0 * k) - q / (6.0 * g);
    let free = q * q / (6.0 * g) + q * q / (18.0 * k) + 0.5 * HARDENING * plastic * plastic;
    // Differentiate free energy at FIXED total F/history, then insert the
    // free-sided stress solution. The temperature derivative is not taken
    // along the force-controlled equilibrium trajectory.
    let entropy = (EXPANSION + EXPANSION_SLOPE * delta_t) * q
        - G_SLOPE * q * q / (6.0 * g * g)
        - K_SLOPE * q * q / (18.0 * k * k);
    let old_g = MU + G_SLOPE * (old.temperature - 300.0);
    let old_k = BULK + K_SLOPE * (old.temperature - 300.0);
    let coefficient_change = (g - old_g) * old.q * old.q / (6.0 * old_g * old_g)
        + 0.5 * (k - old_k) * (old.q / (3.0 * old_k)).powi(2);
    let ts_change = temperature * entropy - old.temperature * old.entropy;
    let thermal_work = q * (old.a - a).exp() * (log_theta - old.log_thermal_stretch).exp_m1();
    let reversible = thermal_work - coefficient_change - ts_change;
    let plastic_heat = YIELD * dp;
    let creep_heat = q * dc;
    let work = -q * (old.a - a).exp_m1();
    let defect =
        work - thermal_work - (free - old.free - coefficient_change) - plastic_heat - creep_heat;
    Reference {
        temperature,
        q,
        a,
        b,
        creep,
        plastic,
        resistance: (old.resistance + dc) / divisor,
        log_thermal_stretch: log_theta,
        free,
        entropy,
        plastic_heat,
        creep_heat,
        coefficient_change,
        ts_change,
        reversible,
        work,
        defect,
    }
}

fn coupled_reference(
    old: Reference,
    dt: f64,
    load: f64,
    external_delta_k: f64,
    combined: bool,
) -> Reference {
    let residual = |temperature: f64| {
        let endpoint = mechanical_reference(old, temperature, dt, load, combined);
        RHO * CP * (temperature - old.temperature - external_delta_k)
            - endpoint.plastic_heat
            - endpoint.creep_heat
            - endpoint.reversible
    };
    let (mut low, mut high) = (
        old.temperature - 5.0,
        old.temperature + external_delta_k + 15.0,
    );
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..100 {
        let temperature = 0.5 * (low + high);
        if residual(temperature) < 0.0 {
            low = temperature;
        } else {
            high = temperature;
        }
    }
    mechanical_reference(old, 0.5 * (low + high), dt, load, combined)
}

#[test]
fn finite_recovering_body_couples_heat_softening_permanent_set_dwell_and_reload() {
    context(|cx| {
        for combined in [false, true] {
            let fixture = recovering_fixture(combined);
            let mut thermal = adiabatic(&fixture);
            let mut solid = FiniteThermoplasticSolid::try_new_with_recovery(
                &thermal,
                &NODES,
                &TETS,
                vec![FiniteThermoplasticCell {
                    cell: CELL,
                    query: fixture.query(),
                }],
                TetAssemblyBudget::default(),
                cx,
            )
            .unwrap();
            let mut observer = IntervalHeat::default();
            let mut old = Reference {
                temperature: 300.0,
                ..Reference::default()
            };
            let mut time = 0.0;
            let force = if combined { 2.2e6 } else { DEAD_PIOLA };
            let prescribed: Vec<_> = (0..12)
                .filter(|dof| ![3, 7, 11].contains(dof))
                .map(|dof| (dof, 0.0))
                .collect();
            for (index, (dt, load, heating_k)) in [
                (0.5, force, 0.0),
                (0.75, force, 15.0),
                (0.25, 0.0, 0.0),
                (2.0, 0.0, 20.0),
                (0.5, force, 0.0),
            ]
            .into_iter()
            .enumerate()
            {
                let expected = coupled_reference(old, dt, load, heating_k, combined);
                let external_j = RHO * V * CP * heating_k;
                observer.watts = external_j / dt;
                let old_delivered = observer.delivered_j;
                let saved_history = solid.mechanics().element_states().to_vec();
                let saved_cells = thermal.cells().clone();
                let forces = [(3, V / L * load)];
                if index == 1 {
                    assert!(
                        solid
                            .advance(
                                &thermal,
                                &observer,
                                time + dt,
                                &prescribed,
                                &forces,
                                [0.0; 3],
                                settings(),
                                controls(),
                                thermal_config(),
                                FiniteThermoplasticConfig {
                                    maximum_iterations: 1,
                                    ..config()
                                },
                                cx,
                            )
                            .is_err()
                    );
                    assert_eq!(solid.mechanics().element_states(), saved_history.as_slice());
                    assert_eq!(thermal.cells(), &saved_cells);
                    close(observer.delivered_j, old_delivered, 0.0);
                }
                let candidate = solid
                    .advance(
                        &thermal,
                        &observer,
                        time + dt,
                        &prescribed,
                        &forces,
                        [0.0; 3],
                        settings(),
                        controls(),
                        thermal_config(),
                        config(),
                        cx,
                    )
                    .unwrap();
                let actual = candidate.solid().mechanics();
                let history = &actual.element_states()[0];
                let report = candidate.report();
                close(actual.temperatures_k()[0], expected.temperature, 3e-7);
                close(
                    candidate.thermal().temperature_k(CELL).unwrap(),
                    expected.temperature,
                    3e-7,
                );
                close(
                    history.native.alpha,
                    if combined {
                        expected.plastic
                    } else {
                        expected.creep
                    },
                    3e-8,
                );
                close(
                    history.native.creep_alpha,
                    if combined { expected.creep } else { 0.0 },
                    3e-8,
                );
                close(history.resistance_strain, expected.resistance, 3e-8);
                let total = expected.plastic + expected.creep;
                for (actual, expected) in history.native.inverse_plastic_cauchy_green.iter().zip([
                    (-2.0 * total).exp(),
                    total.exp(),
                    total.exp(),
                    0.0,
                    0.0,
                    0.0,
                ]) {
                    close(*actual, expected, 4e-8);
                }
                for (axis, strain) in [expected.a, expected.b, expected.b].into_iter().enumerate() {
                    close(
                        actual.displacement_m()[axis + 1][axis],
                        L * strain.exp_m1(),
                        4e-10,
                    );
                    close(
                        actual.element_first_piola_pa()[0][4 * axis],
                        if axis == 0 { load } else { 0.0 },
                        0.5,
                    );
                }
                close(actual.free_energy_j(), V * expected.free, 6e-9);
                close(
                    report.material_energy.plastic_dissipation_j,
                    V * if combined {
                        expected.plastic_heat
                    } else {
                        expected.creep_heat
                    },
                    6e-9,
                );
                close(
                    report.material_energy.creep_dissipation_j,
                    V * if combined { expected.creep_heat } else { 0.0 },
                    6e-9,
                );
                close(
                    report.material_energy.reversible_heat_j,
                    V * expected.reversible,
                    6e-9,
                );
                close(
                    report.material_energy.coefficient_storage_change_j,
                    V * expected.coefficient_change,
                    6e-9,
                );
                close(
                    report.material_energy.thermoelastic_storage_change_j,
                    V * expected.ts_change,
                    6e-9,
                );
                close(
                    report.material_energy.internal_endpoint_work_j,
                    V * expected.work,
                    6e-9,
                );
                close(
                    report.material_energy.endpoint_work_defect_j,
                    V * expected.defect,
                    6e-9,
                );
                close(
                    report.applied_heat_j,
                    V * (expected.plastic_heat + expected.creep_heat + expected.reversible),
                    6e-9,
                );
                close(report.external_heat_j, external_j, 3e-9);
                close(
                    report.thermal_audit.stored_energy_change_j,
                    RHO * V * CP * (expected.temperature - old.temperature),
                    3e-8,
                );
                close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
                close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
                close(
                    candidate.observer().delivered_j,
                    old_delivered + external_j,
                    3e-9,
                );
                assert_eq!(solid.mechanics().element_states(), saved_history.as_slice());
                assert_eq!(thermal.cells(), &saved_cells);
                close(observer.delivered_j, old_delivered, 0.0);
                assert_eq!(
                    actual.families()[0].family().model_card(),
                    &fixture.card.models()[0]
                );
                assert!(std::ptr::eq(
                    actual.families()[0].caloric_source(),
                    &fixture.caloric
                ));
                if index == 3 {
                    close(expected.creep, old.creep, 0.0);
                    close(expected.plastic, old.plastic, 0.0);
                    assert!(expected.resistance < 0.6 * old.resistance);
                    close(report.material_energy.inelastic_dissipation_j(), 0.0, 6e-9);
                    // Thermal substeps do not repeat the mechanical recovery
                    // interval or consume the accepted history more than once.
                    let refined = solid
                        .advance(
                            &thermal,
                            &observer,
                            time + dt,
                            &prescribed,
                            &forces,
                            [0.0; 3],
                            settings(),
                            controls(),
                            StepConfig {
                                max_dt_s: 0.1,
                                ..thermal_config()
                            },
                            config(),
                            cx,
                        )
                        .unwrap();
                    close(
                        refined.solid().mechanics().element_states()[0].resistance_strain,
                        history.resistance_strain,
                        3e-10,
                    );
                    close(
                        refined.thermal().temperature_k(CELL).unwrap(),
                        expected.temperature,
                        3e-7,
                    );
                }
                let (next_thermal, next_solid, next_observer, _) = candidate.into_parts();
                thermal = next_thermal;
                solid = next_solid;
                observer = next_observer;
                old = expected;
                time += dt;
            }
            assert!(solid.mechanics().displacement_m()[1][0] > 0.1 * L);
        }
    });
}
