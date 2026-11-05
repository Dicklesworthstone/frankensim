//! G0 admission/atomicity and G1 independently closed endpoint energy balances.
#[path = "tests/solve.rs"]
mod solve;

#[path = "tests/timeline.rs"]
mod timeline;

#[path = "tests/network.rs"]
mod network;

#[path = "tests/adaptive.rs"]
mod adaptive;

#[path = "tests/activation.rs"]
mod activation;

#[path = "tests/films.rs"]
mod films;

#[path = "tests/resolved_timeline.rs"]
mod resolved_timeline;

use super::*;
use crate::enthalpy::backward_euler::{
    BoundaryExchanges, BoundaryStepCoefficients, PrescribedBathHeat, StepCoefficients,
    TransportStepCoefficients,
};
use crate::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, Event, HeatTransfer, MaterialId,
    Reservoir, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};
use crate::evaporation::{EvaporationError, EvaporationProperties, MOLAR_GAS_CONSTANT};
use std::collections::BTreeMap;

const CP: f64 = 4000.0;
const LATENT: f64 = 2.4e6;

struct Caloric {
    datum: f64,
    fixed_mass: Option<f64>,
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + (h - self.datum) / CP;
        if (280.0..=330.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / CP
    }
    fn validate_cell(&self, cell: &DepositedCell) -> Result<(), DepositionError> {
        if self.fixed_mass.is_some_and(|mass| cell.mass_kg != mass) {
            Err(DepositionError::InvalidInput(
                "original fixed inventory changed",
            ))
        } else {
            Ok(())
        }
    }
}
static CALORIC: Caloric = Caloric {
    datum: 0.0,
    fixed_mass: None,
};

// Analytic interface fixture, not a claimed empirical liquid correlation.
struct Interface {
    pressure_slope: f64,
    latent_slope: f64,
}
impl EvaporationLaw for Interface {
    fn temperature_range_k(&self) -> [f64; 2] {
        [280.0, 330.0]
    }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        Ok(EvaporationProperties {
            saturation_pressure_pa: 4000.0 + self.pressure_slope * (t - 300.0),
            saturation_pressure_derivative_pa_k: self.pressure_slope,
            latent_heat_j_kg: LATENT + self.latent_slope * (t - 300.0),
            liquid_density_kg_m3: 1000.0,
            molar_mass_kg_mol: 0.018,
        })
    }
}
const CONSTANT: Interface = Interface {
    pressure_slope: 0.0,
    latent_slope: 0.0,
};
const VARYING: Interface = Interface {
    pressure_slope: 100.0,
    latent_slope: -2000.0,
};

fn film(vapor_partial_pressure_pa: f64) -> GasFilm {
    GasFilm {
        total_pressure_pa: 100_000.0,
        vapor_partial_pressure_pa,
        temperature_k: 300.0,
        mass_transfer_coefficient_m_s: 0.01,
    }
}
fn transfer(law: &Interface, t: f64, film: GasFilm, area: f64, dt: f64) -> f64 {
    let psat = 4000.0 + law.pressure_slope * (t - 300.0);
    film.mass_transfer_coefficient_m_s * 0.018 * film.total_pressure_pa
        / (MOLAR_GAS_CONSTANT * film.temperature_k)
        * ((film.total_pressure_pa - film.vapor_partial_pressure_pa)
            / (film.total_pressure_pa - psat))
            .ln()
        * area
        * dt
}
fn range(datum: f64) -> [f64; 2] {
    [datum - 20.0 * CP, datum + 30.0 * CP]
}
fn config() -> StepConfig {
    StepConfig {
        max_dt_s: 1000.0,
        energy_tolerance_j: 1e-6,
    }
}
fn cell(id: u64, mass: f64, t: f64, datum: f64, source_w: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: id,
        material_id: MaterialId(1),
        mass_kg: mass,
        specific_h_j_kg: datum + CP * (t - 300.0),
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64, 0.1, 0.0],
            volume_m3: mass / 1000.0,
            surface_area_m2: 0.01,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w,
        },
    }
}
fn simulation<'m>(
    law: &'m dyn ThermalLaw,
    cells: Vec<DepositedCell>,
    link: Option<f64>,
    later: Vec<TimedEvent>,
) -> DepositionSimulation<'m> {
    let mut events: Vec<_> = cells
        .into_iter()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    if let Some(conductance_w_k) = link {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 9,
                a: CellId(1),
                b: CellId(2),
                area_m2: 0.002,
                conductance_w_k,
            }),
        });
    }
    events.extend(later);
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 12,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0
}
fn coefficients<T>(payload: T) -> BoundaryExchangeStepCoefficients<T> {
    TransportStepCoefficients::from(BoundaryStepCoefficients::from(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload,
    }))
    .into()
}
fn empty(
    _: &DepositionSimulation<'_>,
    _: LiquidSurfaceState,
) -> Result<BoundaryExchangeStepCoefficients<()>, ThermalEvaporationError> {
    Ok(coefficients(()))
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn two_cell_backward_euler_retains_conservative_contacts_and_original_baths() {
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
    // Independent 2x2 backward-Euler solution in temperature increments.
    // At constant saturation and L, the signed phase heat is a constant sink.
    let dm = transfer(&CONSTANT, 300.0, film(0.0), 0.006, dt);
    let latent_w = dm * LATENT / dt;
    let diagonal = [400.0 + dt * 2.6, 800.0 + dt * 2.16];
    let rhs = [
        dt * (10.0 + 0.3 * 10.0 + 0.3 * 20.0 + 2.0 * 10.0 - latent_w),
        dt * (-2.0 + 0.16 * -15.0 + 2.0 * -10.0),
    ];
    let determinant = diagonal[0] * diagonal[1] - (dt * 2.0).powi(2);
    let delta = [
        (rhs[0] * diagonal[1] + dt * 2.0 * rhs[1]) / determinant,
        (diagonal[0] * rhs[1] + dt * 2.0 * rhs[0]) / determinant,
    ];
    let expected_t = [300.0 + delta[0], 310.0 + delta[1]];
    assert!(dt > 8.0 * 0.8 * 400.0 / 2.6);
    let step = EvaporatingBackwardEulerStep::new(
        &old,
        CellId(1),
        &CONSTANT,
        film(0.0),
        dt,
        range(0.0),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(step.cell_ids(), &[CellId(1), CellId(2)]);
    assert_eq!(step.duration_s(), dt);
    assert_eq!(step.target_s(), dt);
    let trial = step
        .evaluate(&[400.0 * delta[0], 800.0 * delta[1]], empty, || Ok(()))
        .unwrap();
    assert!(trial.maximum_cell_residual_j() < 1e-7);
    let (next, (), audit) = trial.accept(1e-7, 1e-7, 1e-14, || Ok(())).unwrap();
    for (id, t) in [(CellId(1), expected_t[0]), (CellId(2), expected_t[1])] {
        close(next.temperature_k(id).unwrap(), t, 2e-11);
        assert_eq!(next.cells()[&id].boundary, old.cells()[&id].boundary);
        assert_eq!(next.cells()[&id].geometry, old.cells()[&id].geometry);
    }
    close(audit.surface.transferred_mass_kg, dm, 1e-16);
    let external = dt
        * (10.0 - 2.0
            + 0.3 * (310.0 - expected_t[0])
            + 0.3 * (320.0 - expected_t[0])
            + 0.16 * (295.0 - expected_t[1]));
    close(audit.external_heat_j, external, 1e-7);
    close(
        next.energy_j() - old.energy_j() + audit.surface.vapor_energy_j,
        external,
        1e-7,
    );
    close(next.mass_kg() + dm, old.mass_kg(), 1e-16);
    close(
        audit.stored_energy_change_j,
        next.energy_j() - old.energy_j(),
        1e-7,
    );
    assert_eq!(next.links(), old.links());
    assert_eq!(next.time_s(), dt);
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn variable_latent_heat_and_enthalpy_datum_share_one_open_system_balance() {
    let dt = 1000.0;
    let t = 315.0;
    let dm = transfer(&VARYING, t, film(0.0), 0.01, dt);
    let latent = LATENT - 2000.0 * (t - 300.0);
    let supplied = 0.1 * CP * (t - 300.0) + dm * latent;
    assert!(300.0 + supplied / (0.1 * CP) > 330.0);
    assert!(dm * (LATENT - latent) > 100.0);
    let mut receipts = Vec::new();
    for datum in [0.0, 1e9] {
        let caloric = Caloric {
            datum,
            fixed_mass: None,
        };
        let old = simulation(
            &caloric,
            vec![cell(1, 0.1, 300.0, datum, supplied / dt)],
            None,
            vec![],
        );
        let step = EvaporatingBackwardEulerStep::new(
            &old,
            CellId(1),
            &VARYING,
            film(0.0),
            dt,
            range(datum),
            || Ok(()),
        )
        .unwrap();
        let trial = step.evaluate(&[0.1 * CP * 15.0], empty, || Ok(())).unwrap();
        close(trial.audit().surface.latent_energy_j, dm * latent, 1e-7);
        let (next, (), audit) = trial.accept(1e-7, 1e-7, 1e-14, || Ok(())).unwrap();
        close(next.temperature_k(CellId(1)).unwrap(), t, 1e-12);
        close(audit.external_heat_j, supplied, 1e-10);
        close(
            next.energy_j() - old.energy_j() + audit.surface.vapor_energy_j,
            supplied,
            1e-6,
        );
        close(audit.energy_residual_j, 0.0, 1e-7);
        receipts.push(audit);
    }
    assert_eq!(
        receipts[0].surface.liquid.mass_kg,
        receipts[1].surface.liquid.mass_kg
    );
    close(
        receipts[1].surface.vapor_energy_j - receipts[0].surface.vapor_energy_j,
        receipts[0].surface.transferred_mass_kg * 1e9,
        1e-8,
    );
    close(
        receipts[1].stored_energy_change_j - receipts[0].stored_energy_change_j,
        -receipts[0].surface.transferred_mass_kg * 1e9,
        1e-8,
    );
}

#[test]
fn endpoint_callback_uses_new_inventory_and_original_boundary_exchange_ingress() {
    let dt = 1000.0;
    let target = [305.0, 302.0];
    let mut a = cell(1, 0.1, 300.0, 0.0, 3.0);
    a.boundary.reservoir = Some(Reservoir {
        area_m2: 0.002,
        conductance_w_k: 7.0,
        temperature_k: 320.0,
    });
    let mut b = cell(2, 0.2, 300.0, 0.0, 0.0);
    b.boundary.ambient_k = 290.0;
    let old = simulation(&CALORIC, vec![a, b], Some(0.25), vec![]);
    let dm = transfer(&VARYING, target[0], film(0.0), 0.006, dt);
    let mass = 0.1 - dm;
    let latent = LATENT - 2000.0 * 5.0;
    let baths = |t: f64| {
        [
            PrescribedBathHeat {
                outward_power_w: 5.670_374_419e-8 * 0.006 * (t.powi(4) - 280.0_f64.powi(4)),
                temperature_k: 280.0,
            },
            PrescribedBathHeat {
                outward_power_w: 0.02 * (t - 320.0),
                temperature_k: 320.0,
            },
        ]
    };
    let outward = baths(target[0])
        .iter()
        .map(|bath| bath.outward_power_w)
        .sum::<f64>();
    let liquid_heat = 0.1 * CP * 5.0 / dt + dm * latent / dt;
    let base_source = liquid_heat - (3.0 - outward + 0.4 * 15.0 - 1.5 * 3.0 - 0.75) - 20.0 * mass;
    let other_source = 0.2 * CP * 2.0 / dt - (0.1 * -12.0 + 1.5 * 3.0 + 0.75);
    let step = EvaporatingBackwardEulerStep::new(
        &old,
        CellId(1),
        &VARYING,
        film(0.0),
        dt,
        range(0.0),
        || Ok(()),
    )
    .unwrap();
    struct Receipt {
        proposed_mass: f64,
        original_mass: f64,
    }
    let trial = step
        .evaluate(
            &[2000.0, 1600.0],
            |thermal, liquid| {
                close(thermal.temperature_k(CellId(1))?, target[0], 1e-12);
                assert_eq!(thermal.time_s(), dt);
                let mut result = coefficients(Receipt {
                    proposed_mass: liquid.mass_kg,
                    original_mass: thermal.cells()[&CellId(1)].mass_kg,
                });
                result.transport.heat.bulk.sources_w = BTreeMap::from([
                    (CellId(1), base_source + 20.0 * liquid.mass_kg),
                    (CellId(2), other_source),
                ]);
                result.transport.heat.bulk.conductances_w_k.insert(9, 1.5);
                result
                    .transport
                    .heat
                    .boundaries
                    .ambient_w_k
                    .insert(CellId(2), 0.1);
                result
                    .transport
                    .heat
                    .boundaries
                    .reservoir_w_k
                    .insert(CellId(1), 0.4);
                result.transport.transfers.push(HeatTransfer {
                    a: CellId(1),
                    b: CellId(2),
                    power_w: 0.75,
                });
                result.exchanges = BoundaryExchanges {
                    ambient: BTreeMap::from([(
                        CellId(1),
                        baths(thermal.temperature_k(CellId(1))?),
                    )]),
                    reservoir: BTreeMap::new(),
                };
                Ok::<_, ThermalEvaporationError>(result)
            },
            || Ok(()),
        )
        .unwrap();
    assert_eq!(trial.payload().original_mass, 0.1);
    close(trial.payload().proposed_mass, mass, 1e-16);
    let (next, payload, audit) = trial.accept(1e-7, 1e-7, 1e-14, || Ok(())).unwrap();
    assert!(payload.proposed_mass < payload.original_mass);
    close(audit.surface.supplied_heat_j, liquid_heat * dt, 1e-7);
    let external =
        dt * (3.0 + base_source + 20.0 * mass + other_source - outward + 0.4 * 15.0 + 0.1 * -12.0);
    close(audit.external_heat_j, external, 1e-7);
    close(
        next.energy_j() + audit.surface.vapor_energy_j,
        external,
        1e-7,
    );
    assert_eq!(next.links(), old.links());
    for id in step.cell_ids() {
        assert_eq!(next.cells()[id].boundary, old.cells()[id].boundary);
    }
}

#[test]
fn heating_and_cooling_reverse_the_phase_flux_at_the_actual_endpoint() {
    for (t, pressure) in [(306.0, 4500.0), (294.0, 3500.0)] {
        let dt = 2000.0;
        let film = film(pressure);
        let dm = transfer(&VARYING, t, film, 0.01, dt);
        assert!(dm * transfer(&VARYING, 300.0, film, 0.01, dt) < 0.0);
        let latent = LATENT - 2000.0 * (t - 300.0);
        let q = 0.1 * CP * (t - 300.0) + dm * latent;
        let old = simulation(
            &CALORIC,
            vec![cell(1, 0.1, 300.0, 0.0, q / dt)],
            None,
            vec![],
        );
        let step = EvaporatingBackwardEulerStep::new(
            &old,
            CellId(1),
            &VARYING,
            film,
            dt,
            range(0.0),
            || Ok(()),
        )
        .unwrap();
        let trial = step
            .evaluate(&[0.1 * CP * (t - 300.0)], empty, || Ok(()))
            .unwrap();
        let (next, (), audit) = trial.accept(1e-7, 1e-7, 1e-14, || Ok(())).unwrap();
        close(next.cells()[&CellId(1)].mass_kg, 0.1 - dm, 1e-16);
        close(audit.surface.transferred_mass_kg, dm, 1e-16);
        close(next.energy_j() + audit.surface.vapor_energy_j, q, 1e-7);
        assert_eq!(
            audit.surface.transferred_mass_kg.is_sign_positive(),
            t > 300.0
        );
    }
}

#[test]
fn original_inventory_source_and_positive_mass_admission_are_retained() {
    let fixed = Caloric {
        datum: 0.0,
        fixed_mass: Some(0.1),
    };
    let old = simulation(&fixed, vec![cell(1, 0.1, 300.0, 0.0, 0.0)], None, vec![]);
    let step = EvaporatingBackwardEulerStep::new(
        &old,
        CellId(1),
        &VARYING,
        film(0.0),
        1.0,
        range(0.0),
        || Ok(()),
    )
    .unwrap();
    assert!(matches!(
        step.evaluate(&[0.0], empty, || Ok(())),
        Err(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput("original fixed inventory changed")
        ))
    ));
    let old = simulation(&CALORIC, vec![cell(1, 0.1, 300.0, 0.0, 0.0)], None, vec![]);
    for dt in [1e6, 1e-20] {
        let step = EvaporatingBackwardEulerStep::new(
            &old,
            CellId(1),
            &VARYING,
            film(0.0),
            dt,
            range(0.0),
            || Ok(()),
        )
        .unwrap();
        match step.evaluate(&[0.0], empty, || Ok(())) {
            Err(ThermalEvaporationError::Surface(EvaporationError::InventoryLimit {
                maximum_duration_s,
            })) if dt == 1e6 => assert!(maximum_duration_s > 0.0 && maximum_duration_s < dt),
            Err(ThermalEvaporationError::Surface(EvaporationError::Unrepresentable(_)))
                if dt == 1e-20 => {}
            _ => panic!("wrong inventory/representability admission"),
        }
    }
    assert!(matches!(
        EvaporatingBackwardEulerStep::new(
            &old,
            CellId(1),
            &VARYING,
            GasFilm {
                total_pressure_pa: 4000.0,
                ..film(0.0)
            },
            1.0,
            range(0.0),
            || Ok(()),
        ),
        Err(ThermalEvaporationError::Surface(
            EvaporationError::Boiling { .. }
        ))
    ));
    let step = EvaporatingBackwardEulerStep::new(
        &old,
        CellId(1),
        &VARYING,
        film(0.0),
        1.0,
        [0.0, CP],
        || Ok(()),
    )
    .unwrap();
    assert!(matches!(
        step.evaluate(&[0.1 * CP * 2.0], empty, || Ok(())),
        Err(ThermalEvaporationError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
    assert!(matches!(
        step.evaluate(
            &[0.0],
            |_, _| {
                Err::<BoundaryExchangeStepCoefficients<()>, _>(ThermalEvaporationError::Surface(
                    EvaporationError::Invalid("original source callback refused"),
                ))
            },
            || Ok(())
        ),
        Err(ThermalEvaporationError::Surface(EvaporationError::Invalid(
            "original source callback refused"
        )))
    ));
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.cells()[&CellId(1)].mass_kg, 0.1);
    assert_eq!(old.cells()[&CellId(1)].specific_h_j_kg, 0.0);
}

#[test]
fn local_mass_and_final_cancellation_gates_cannot_be_bypassed_by_total_balance() {
    let old = simulation(
        &CALORIC,
        vec![
            cell(1, 0.1, 300.0, 0.0, 0.0),
            cell(2, 0.1, 300.0, 0.0, 0.0),
            cell(3, 0.1, 300.0, 0.0, 0.0),
        ],
        None,
        vec![],
    );
    let sealed = GasFilm {
        mass_transfer_coefficient_m_s: 0.0,
        ..film(0.0)
    };
    let step = EvaporatingBackwardEulerStep::new(
        &old,
        CellId(1),
        &VARYING,
        sealed,
        1.0,
        range(0.0),
        || Ok(()),
    )
    .unwrap();
    let trial = step.evaluate(&[0.0, 1.0, -1.0], empty, || Ok(())).unwrap();
    assert_eq!(trial.audit().energy_residual_j, 0.0);
    assert_eq!(trial.audit().surface.energy_residual_j, 0.0);
    assert_eq!(trial.maximum_cell_residual_j(), 1.0);
    assert!(matches!(
        trial.accept(1e-9, 1e-9, 1e-14, || Ok(())),
        Err(ThermalEvaporationError::Thermal(
            DepositionError::EnergyBalance(_)
        ))
    ));
    let mut polls = 0;
    step.evaluate(&[0.0; 3], empty, || Ok(()))
        .unwrap()
        .accept(1e-9, 1e-9, 1e-14, || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    for cancel_at in [1, polls] {
        let trial = step.evaluate(&[0.0; 3], empty, || Ok(())).unwrap();
        let mut count = 0;
        assert!(matches!(
            trial.accept(1e-9, 1e-9, 1e-14, || {
                count += 1;
                if count == cancel_at {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(ThermalEvaporationError::Thermal(DepositionError::Cancelled))
        ));
    }
    let q = transfer(&CONSTANT, 300.0, film(0.0), 0.01, 1.0) * LATENT;
    let old = simulation(&CALORIC, vec![cell(1, 0.1, 300.0, 0.0, q)], None, vec![]);
    let step = EvaporatingBackwardEulerStep::new(
        &old,
        CellId(1),
        &CONSTANT,
        film(0.0),
        1.0,
        range(0.0),
        || Ok(()),
    )
    .unwrap();
    let trial = step.evaluate(&[0.0], empty, || Ok(())).unwrap();
    assert!(trial.audit().surface.mass_rate_residual_kg != 0.0);
    assert!(matches!(
        trial.accept(1e-7, 1e-7, 0.0, || Ok(())),
        Err(ThermalEvaporationError::MassBalance { .. })
    ));
    assert_eq!(old.cells()[&CellId(1)].mass_kg, 0.1);
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn endpoint_events_remain_pending_until_the_coupled_interval_is_accepted() {
    let power = transfer(&CONSTANT, 300.0, film(0.0), 0.01, 1.0) * LATENT;
    let initial = cell(1, 0.1, 300.0, 0.0, power);
    let mut covered = initial.boundary;
    covered.source_w = 0.0;
    covered.reservoir = Some(Reservoir {
        area_m2: 0.01,
        temperature_k: 300.0,
        conductance_w_k: 0.0,
    });
    let old = simulation(
        &CALORIC,
        vec![initial],
        None,
        vec![TimedEvent {
            time_s: 10.0,
            event: Event::Boundary(CellId(1), covered),
        }],
    );
    assert!(matches!(
        EvaporatingBackwardEulerStep::new(
            &old,
            CellId(1),
            &CONSTANT,
            film(0.0),
            11.0,
            range(0.0),
            || Ok(()),
        ),
        Err(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput(_)
        ))
    ));
    let step = EvaporatingBackwardEulerStep::new(
        &old,
        CellId(1),
        &CONSTANT,
        film(0.0),
        10.0,
        range(0.0),
        || Ok(()),
    )
    .unwrap();
    let (left, (), audit) = step
        .evaluate(&[0.0], empty, || Ok(()))
        .unwrap()
        .accept(1e-7, 1e-7, 1e-14, || Ok(()))
        .unwrap();
    assert_eq!(left.time_s(), 10.0);
    assert_eq!(
        left.cells()[&CellId(1)].boundary,
        old.cells()[&CellId(1)].boundary
    );
    close(
        audit.surface.transferred_mass_kg,
        transfer(&CONSTANT, 300.0, film(0.0), 0.01, 10.0),
        1e-16,
    );
    assert!(matches!(
        EvaporatingBackwardEulerStep::new(
            &left,
            CellId(1),
            &CONSTANT,
            film(0.0),
            20.0,
            range(0.0),
            || Ok(()),
        ),
        Err(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput(_)
        ))
    ));
    let (right, _) = left.advance_to(10.0, config(), || Ok(())).unwrap();
    assert_eq!(right.cells()[&CellId(1)].boundary, covered);
    assert_eq!(right.exposed_area_m2(CellId(1)).unwrap(), 0.0);
    assert_eq!(
        right.cells()[&CellId(1)].mass_kg,
        left.cells()[&CellId(1)].mass_kg
    );
}
