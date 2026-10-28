//! Read-only half-space oracle versus the original finite-volume owner.
//! Synthetic constant materials; no moving-source validity is inferred.
use fs_couple::friction_heat::flash::{FlashHeatHistory, FlashOracleError};
use fs_couple::friction_heat::{
    FrictionHeatInput, FrictionHeatReport, FrictionHeatRoute, FrictionHeatStep,
    FrictionHeatTolerances,
};
use fs_la::prepared_lu::LuWorkspace;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::{
    FlashTemperatureEstimate, FlashTemperatureInput, HeatPartition, InputAuthority,
    InterfaceMedium, InterfaceSystemRef, SurfaceThermalProperties, flash_temperature_candidate,
};
use std::cell::Cell;
use std::collections::BTreeMap;
struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h;
        if t.is_finite() && t > 0.0 {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
static CALORIC: Caloric = Caloric;
const CONTACT: u64 = 0;
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t, "{a:.16e} != {b:.16e}");
}
fn interface() -> InterfaceSystemRef {
    InterfaceSystemRef::new(
        "synthetic/halfspace-a->b",
        "synthetic/initial-uniform",
        "synthetic/constant-caloric-conductivity",
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .unwrap()
}
fn properties(k: [f64; 2]) -> [SurfaceThermalProperties; 2] {
    // rho*cp=1 J/(m3 K), hence alpha=k in this synthetic fixture.
    k.map(|v| SurfaceThermalProperties::new(v, v).unwrap())
}
fn network(n: usize, k: [f64; 2], initial: [f64; 2]) -> DepositionSimulation<'static> {
    let dx = 4.0 / n as f64;
    let mut events = Vec::new();
    for side in 0..2 {
        for i in 0..n {
            let id = CellId((side * n + i + 1) as u64);
            let (start, end) = if side == 0 {
                (-((i + 1) as f64) * dx, -(i as f64) * dx)
            } else {
                (i as f64 * dx, (i + 1) as f64 * dx)
            };
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id,
                    segment_id: id.0,
                    body_id: side as u64,
                    material_id: MaterialId(0),
                    mass_kg: dx,
                    specific_h_j_kg: initial[side] - 300.0,
                    geometry: CellGeometry {
                        start_m: [start, 0.0, 0.0],
                        end_m: [end, 0.0, 0.0],
                        volume_m3: dx,
                        surface_area_m2: 2.0 + 4.0 * dx,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            });
        }
        for i in 0..n - 1 {
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: (side * n + i + 1) as u64,
                    a: CellId((side * n + i + 1) as u64),
                    b: CellId((side * n + i + 2) as u64),
                    area_m2: 1.0,
                    conductance_w_k: k[side] / dx,
                }),
            });
        }
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: CONTACT,
            a: CellId(1),
            b: CellId((n + 1) as u64),
            area_m2: 1.0,
            conductance_w_k: series_conductance_w_k(1.0, [0.5 * dx; 2], k).unwrap(),
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), &CALORIC as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2 * n,
            max_links: 2 * n,
            max_events: 4 * n + 2,
            max_substeps: 4096,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn factor(thermal: &DepositionSimulation<'_>, dt: f64, k: [f64; 2]) -> LuWorkspace {
    let count = thermal.cells().len();
    let dimension = count + 1;
    let mut matrix = vec![0.0; dimension * dimension];
    for (i, cell) in thermal.cells().values().enumerate() {
        matrix[i * dimension + i] = cell.mass_kg / dt;
    }
    for link in thermal.links().values().filter(|link| link.id != CONTACT) {
        let a = link.a.0 as usize - 1;
        let b = link.b.0 as usize - 1;
        let g = link.conductance_w_k;
        matrix[a * dimension + a] += g;
        matrix[b * dimension + b] += g;
        matrix[a * dimension + b] -= g;
        matrix[b * dimension + a] -= g;
    }
    for (side, i) in [0, count / 2].into_iter().enumerate() {
        let g = k[side] * count as f64 / 4.0;
        matrix[i * dimension + i] += g;
        matrix[i * dimension + count] -= g;
        matrix[count * dimension + i] -= g;
        matrix[count * dimension + count] += g;
    }
    let mut lu = LuWorkspace::new(dimension).unwrap();
    lu.factor(&matrix).unwrap();
    lu
}
fn advance(
    thermal: &DepositionSimulation<'static>,
    route: &FrictionHeatRoute,
    lu: &mut LuWorkspace,
    dt: f64,
    power: f64,
    k: [f64; 2],
) -> (DepositionSimulation<'static>, FrictionHeatReport) {
    let mut rhs: Vec<_> = thermal
        .cells()
        .values()
        .map(|cell| cell.mass_kg / dt * cell.specific_h_j_kg)
        .collect();
    rhs.push(power);
    let mut temperatures = vec![0.0; rhs.len()];
    lu.solve_factored_into(&rhs, &mut temperatures).unwrap();
    let energy: Vec<_> = thermal
        .cells()
        .values()
        .zip(&temperatures)
        .map(|(cell, &t)| cell.mass_kg * (t - cell.specific_h_j_kg))
        .collect();
    let step = FrictionHeatStep::new(thermal, thermal.time_s() + dt, route, || Ok(())).unwrap();
    let trial = step
        .evaluate(
            &energy,
            300.0 + temperatures.last().unwrap(),
            |_, _| {
                Ok::<_, DepositionError>(FrictionHeatInput {
                    dissipated_j: power * dt,
                    conductivity_w_m_k: k,
                    transport: StepCoefficients {
                        sources_w: BTreeMap::new(),
                        conductances_w_k: BTreeMap::new(),
                        payload: (),
                    },
                })
            },
            || Ok(()),
        )
        .unwrap();
    let (next, (), report) = trial
        .accept(
            FrictionHeatTolerances {
                cell_j: 1e-10,
                thermal_j: 1e-10,
                interface_j: 1e-10,
                total_j: 1e-10,
            },
            || Ok(()),
        )
        .unwrap();
    (next, report)
}
#[test]
fn g0_explicit_duration_uses_original_flash_response_without_fake_motion() {
    let p = properties([1.0, 4.0]);
    let partition = HeatPartition::semi_infinite_surfaces(p[0], p[1]).unwrap();
    let old =
        FlashTemperatureInput::new(3.0, partition, 1.0, 2.0, 8.0, Some(p[0]), Some(p[1])).unwrap();
    let direct =
        FlashTemperatureInput::for_duration(3.0, partition, 1.0, 0.25, Some(p[0]), Some(p[1]))
            .unwrap();
    let old = flash_temperature_candidate(&interface(), old).unwrap();
    let direct = flash_temperature_candidate(&interface(), direct).unwrap();
    assert_eq!(old, direct);
    let FlashTemperatureEstimate::Candidate(candidate) = direct else {
        panic!("known sources");
    };
    close(
        candidate.surface_a_rise_k(),
        1.0 / std::f64::consts::PI.sqrt(),
        1e-14,
    );
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            FlashTemperatureInput::for_duration(3.0, partition, 1.0, bad, Some(p[0]), Some(p[1]))
                .is_err()
        );
    }
}
#[test]
fn g1_original_flash_oracle_matches_resolved_network_on_three_spatial_rungs() {
    let k = [1.0, 4.0];
    let dt = 0.25 / 512.0;
    let mut errors = Vec::new();
    for n in [16, 32, 64] {
        let mut thermal = network(n, k, [300.0; 2]);
        let route = FrictionHeatRoute::from_link(&thermal, CONTACT, k).unwrap();
        let mut history = FlashHeatHistory::new(
            &thermal,
            &route,
            interface(),
            properties(k),
            2 * n,
            512,
            1e-10,
            || Ok(()),
        )
        .unwrap();
        let mut lu = factor(&thermal, dt, k);
        let mut last = None;
        for _ in 0..512 {
            let (next, report) = advance(&thermal, &route, &mut lu, dt, 3.0, k);
            last = Some(
                history
                    .observe(&thermal, &next, &report, || Ok(()))
                    .unwrap(),
            );
            thermal = next;
        }
        let comparison = last.unwrap();
        let error = comparison
            .relative_rise_error
            .into_iter()
            .map(Option::unwrap)
            .fold(0.0, f64::max);
        println!("flash grid={n} maximum_relative_rise_error={error:.9e}");
        errors.push(error);
        close(thermal.energy_j(), 0.75, 1e-10);
        assert_eq!(history.time_s(), thermal.time_s());
    }
    assert!(errors[1] < errors[0] && errors[2] < errors[1]);
    assert!(
        errors[2] < 0.02,
        "finest resolved network must agree within 2%: {errors:?}"
    );
}
#[test]
fn g1_piecewise_flux_and_cooldown_retain_history_and_match_duhamel() {
    let n = 64;
    let k = [1.0; 2];
    let dt = 0.25 / 256.0;
    let mut thermal = network(n, k, [300.0; 2]);
    let route = FrictionHeatRoute::from_link(&thermal, CONTACT, k).unwrap();
    let mut history = FlashHeatHistory::new(
        &thermal,
        &route,
        interface(),
        properties(k),
        2 * n,
        256,
        1e-10,
        || Ok(()),
    )
    .unwrap();
    let mut lu = factor(&thermal, dt, k);
    let mut peak = 0.0;
    let mut last = None;
    for i in 0..256 {
        let power = if i < 64 {
            3.0
        } else if i < 128 {
            1.0
        } else {
            0.0
        };
        let (next, report) = advance(&thermal, &route, &mut lu, dt, power, k);
        let comparison = history
            .observe(&thermal, &next, &report, || Ok(()))
            .unwrap();
        if i == 127 {
            peak = comparison.oracle_surface_temperature_k[0];
        }
        last = Some(comparison);
        thermal = next;
    }
    let last = last.unwrap();
    let t = 0.25_f64;
    let expected = (3.0 * t.sqrt() - 2.0 * (t - 0.0625).sqrt() - (t - 0.125).sqrt())
        / std::f64::consts::PI.sqrt();
    for side in 0..2 {
        close(
            last.oracle_surface_temperature_k[side] - 300.0,
            expected,
            1e-8,
        );
        assert!(last.relative_rise_error[side].unwrap() < 0.02);
    }
    assert!(last.oracle_surface_temperature_k[0] < peak);
    assert!(history.retained_pulses() > 0);
    close(thermal.energy_j(), 0.25, 1e-10);
}
#[test]
fn g0_signed_exchange_budget_and_cancel_preserve_observer_and_thermal_state() {
    let k = [1.0; 2];
    let dt = 0.01;
    let thermal = network(8, k, [350.0, 300.0]);
    let route = FrictionHeatRoute::from_link(&thermal, CONTACT, k).unwrap();
    let mut history = FlashHeatHistory::new(
        &thermal,
        &route,
        interface(),
        properties(k),
        16,
        1,
        1e-9,
        || Ok(()),
    )
    .unwrap();
    let mut lu = factor(&thermal, dt, k);
    let (next, report) = advance(&thermal, &route, &mut lu, dt, 0.0, k);
    assert!(report.body_heat_j[0] < 0.0 && report.body_heat_j[1] > 0.0);
    let polls = Cell::new(0);
    let cancelled = history.observe(&thermal, &next, &report, || {
        polls.set(polls.get() + 1);
        if polls.get() == 3 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        cancelled,
        Err(FlashOracleError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(history.time_s(), 0.0);
    assert_eq!(history.retained_pulses(), 0);
    let mut replay = history.clone();
    let first = history
        .observe(&thermal, &next, &report, || Ok(()))
        .unwrap();
    let repeated = replay.observe(&thermal, &next, &report, || Ok(())).unwrap();
    assert_eq!(
        first.oracle_surface_temperature_k,
        repeated.oracle_surface_temperature_k
    );
    assert!(
        first.oracle_surface_temperature_k[0] < 350.0
            && first.oracle_surface_temperature_k[1] > 300.0
    );
    let (end, second) = advance(&next, &route, &mut lu, dt, 3.0, k);
    assert!(matches!(
        history.observe(&next, &end, &second, || Ok(())),
        Err(FlashOracleError::HistoryBudget { maximum: 1 })
    ));
    assert_eq!(history.time_s(), dt);
    assert_eq!(history.retained_pulses(), 1);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(next.time_s(), dt);
    let mut wrong = second.clone();
    wrong.conductances_w_k[0] *= 2.0;
    assert!(replay.observe(&next, &end, &wrong, || Ok(())).is_err());
    assert_eq!(replay.time_s(), dt);
    assert!(
        FlashHeatHistory::new(
            &thermal,
            &route,
            interface(),
            properties([2.0, 1.0]),
            16,
            8,
            1e-9,
            || Ok(())
        )
        .is_err()
    );
}
