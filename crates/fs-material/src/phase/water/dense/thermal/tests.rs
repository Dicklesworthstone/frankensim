use super::*;
use crate::phase::water::dense::DenseWaterControls;
use fs_qty::{Density, Temperature};
use fs_thermal::enthalpy::advection::{Advected, MassBoundary};
use fs_thermal::enthalpy::{Boundary, CellGeometry, Reservoir, ThermalLink};
use fs_thermochem::water::if97::region3::increments::DenseIsochoricIntervalSolve;

fn law() -> DenseIsochoricThermalLaw {
    DenseIsochoricThermalLaw::new(
        DenseIsochore::new(
            Density::new(500.0),
            DenseWaterControls {
                minimum_temperature_k: 650.0,
                maximum_temperature_k: 750.0,
                internal_energy_tolerance_j_kg: 1e-6,
                maximum_iterations: 96,
                interval: DenseIsochoricIntervalSolve {
                    maximum_intervals: 65_535,
                },
            },
            || false,
        )
        .unwrap(),
    )
}
fn state(law: &DenseIsochoricThermalLaw, t: f64) -> DenseWaterState {
    law.chart()
        .state_at_temperature(Temperature::new(t), || false)
        .unwrap()
}
fn boundary(bath: Option<(f64, f64)>) -> Boundary {
    Boundary {
        ambient_k: 650.0,
        htc_w_m2_k: 0.0,
        source_w: 0.0,
        reservoir: bath.map(|(t, g)| Reservoir {
            area_m2: 0.01,
            conductance_w_k: g,
            temperature_k: t,
        }),
    }
}
fn birth(
    law: &DenseIsochoricThermalLaw,
    id: u64,
    t: f64,
    mass: f64,
    bath: Option<(f64, f64)>,
) -> TimedEvent {
    TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(id),
            segment_id: id,
            body_id: id,
            material_id: MaterialId(0),
            mass_kg: mass,
            specific_h_j_kg: state(law, t).specific_internal_energy_j_kg(),
            geometry: CellGeometry {
                start_m: [id as f64, 0.0, 0.0],
                end_m: [id as f64, 0.1, 0.0],
                volume_m3: mass / law.chart().density().value(),
                surface_area_m2: 0.02,
            },
            boundary: boundary(bath),
        }),
    }
}
fn budget() -> DepositionBudget {
    DepositionBudget {
        max_cells: 2,
        max_links: 1,
        max_events: 4,
        max_substeps: 8192,
    }
}
fn controls(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 2e-6,
    }
}
fn near(a: f64, b: f64, allowance: f64) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= allowance,
        "{a:.15e} != {b:.15e}; allowance {allowance:.3e}"
    );
}

#[test]
fn source_cv_entropy_and_fixed_volume_binding_cannot_become_advected_enthalpy() {
    let law = law();
    assert_eq!(law.storage_basis(), ThermalStorageBasis::InternalEnergy);
    for t in [650.0, 673.15, 700.0, 750.0] {
        let original = state(&law, t);
        let u = original.specific_internal_energy_j_kg();
        near(law.temperature_k(u).unwrap(), t, 1e-9);
        let tangent = law
            .temperature_internal_energy_derivative_kg_k_j(u, || false)
            .unwrap();
        near(tangent, 1.0 / original.source().cv().value(), 1e-14);
        assert!(tangent <= law.maximum_temperature_slope());
    }
    let u = state(&law, 700.0).specific_internal_energy_j_kg();
    let small_entropy = law
        .specific_entropy_change_j_kg_k(u, u.next_up())
        .unwrap()
        .unwrap();
    assert!(small_entropy > 0.0);
    near(small_entropy / ((u.next_up() - u) / 700.0), 1.0, 1e-10);
    assert_eq!(
        law.specific_entropy_change_j_kg_k(u.next_up(), u)
            .unwrap()
            .unwrap(),
        -small_entropy
    );

    let mut wrong_volume = birth(&law, 0, 700.0, 0.05, None);
    if let Event::Birth(cell) = &mut wrong_volume.event {
        cell.geometry.volume_m3 *= 1.001;
    }
    assert!(matches!(
        DenseWaterThermalNetwork::new(&law, vec![wrong_volume], budget(), || Ok(())),
        Err(DenseWaterThermalError::Thermal(
            DepositionError::InvalidInput(_)
        ))
    ));
    let original = DenseWaterThermalNetwork::new(
        &law,
        vec![birth(&law, 0, 700.0, 0.05, None)],
        budget(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(
        original
            .thermal()
            .material_storage_basis(MaterialId(0))
            .unwrap(),
        ThermalStorageBasis::InternalEnergy
    );
    let open = Advected::new((), vec![], 0)
        .unwrap()
        .with_boundaries(
            vec![
                MassBoundary::Inlet {
                    id: 0,
                    cell: CellId(0),
                    kg_s: 0.001,
                    material_id: MaterialId(0),
                    specific_h_j_kg: u,
                },
                MassBoundary::Outlet {
                    id: 1,
                    cell: CellId(0),
                    kg_s: 0.001,
                },
            ],
            2,
        )
        .unwrap();
    assert!(matches!(
        original
            .thermal()
            .advance_observed(&open, 1.0, controls(0.1), || Ok(())),
        Err(DepositionError::InvalidInput(_))
    ));
    assert_eq!(original.time_s(), 0.0);
}

#[test]
fn bath_transient_converges_to_independent_cv_integral_and_preserves_u_h_pressure_accounts() {
    let law = law();
    let mass = 0.05;
    let conductance = 5.0;
    let bath = 750.0;
    let initial = DenseWaterThermalNetwork::new(
        &law,
        vec![birth(&law, 0, 650.0, mass, Some((bath, conductance)))],
        budget(),
        || Ok(()),
    )
    .unwrap();
    // Independent separation of variables: t=m/G * integral Cv(T)/(Tb-T)dT.
    // This evaluates the original source at temperature, never its u inverse.
    let intervals = 512;
    let width = 50.0 / intervals as f64;
    let mut quadrature = 0.0;
    for i in 0..=intervals {
        let t = 650.0 + width * i as f64;
        let weight = if i == 0 || i == intervals {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        quadrature += weight * state(&law, t).source().cv().value() / (bath - t);
    }
    let target = (mass / conductance) * width * quadrature / 3.0;
    let mut errors = Vec::new();
    for dt in [0.4, 0.2, 0.1] {
        let (next, audit) = initial.advance_to(target, controls(dt), || Ok(())).unwrap();
        let current = next.profile(|| Ok(())).unwrap()[0];
        errors.push((current.material.temperature_k() - 700.0).abs());
        let before = state(&law, 650.0);
        let heat = audit.thermal.external_heat_j;
        near(
            mass * (current.material.specific_internal_energy_j_kg()
                - before.specific_internal_energy_j_kg()),
            heat,
            2e-6,
        );
        near(
            mass * (current.material.specific_enthalpy_j_kg() - before.specific_enthalpy_j_kg())
                - current.volume_m3 * (current.material.pressure_pa() - before.pressure_pa()),
            heat,
            2e-6,
        );
        assert_eq!(audit.boundary_work_j, 0.0);
        assert!(audit.pressure_enthalpy_change_j > 1000.0);
        assert!(current.material.pressure_pa() > before.pressure_pa());
        near(
            audit.entropy.reservoir_entropy_input_j_k,
            heat / bath,
            1e-11,
        );
        assert!(audit.source_entropy_change_j_k - heat / bath > 0.0);
        assert!(audit.entropy.conduction_generation_j_k > 0.0);
        assert!(!audit.entropy.has_unclassified_sources);
        near(
            audit.entropy.thermal_storage_change_j_k - audit.source_entropy_change_j_k,
            audit.entropy_storage_source_defect_j_k,
            1e-14,
        );
    }
    eprintln!("rigid-water bath time={target:.9}s, temperature errors={errors:?} K");
    assert!(errors[1] < errors[0] && errors[2] < errors[1]);
    assert!(errors[2] < 0.4 * errors[0]);
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn sealed_finite_bodies_equilibrate_by_internal_energy_without_spurious_pdv_work() {
    let law = law();
    let masses = [0.02, 0.01];
    let initial = DenseWaterThermalNetwork::new(
        &law,
        vec![
            birth(&law, 0, 650.0, masses[0], None),
            birth(&law, 1, 750.0, masses[1], None),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 0,
                    a: CellId(0),
                    b: CellId(1),
                    area_m2: 0.01,
                    conductance_w_k: 10.0,
                }),
            },
        ],
        budget(),
        || Ok(()),
    )
    .unwrap();
    let before = initial.profile(|| Ok(())).unwrap();
    let target_u = before
        .iter()
        .map(|c| c.mass_kg * c.material.specific_internal_energy_j_kg())
        .sum::<f64>()
        / masses.iter().sum::<f64>();
    let (mut lo, mut hi) = (650.0, 750.0);
    for _ in 0..48 {
        let midpoint = 0.5 * (lo + hi);
        if state(&law, midpoint).specific_internal_energy_j_kg() < target_u {
            lo = midpoint;
        } else {
            hi = midpoint;
        }
    }
    let equilibrium = 0.5 * (lo + hi);
    let (next, audit) = initial.advance_to(100.0, controls(0.1), || Ok(())).unwrap();
    let after = next.profile(|| Ok(())).unwrap();
    for cell in &after {
        near(cell.material.temperature_k(), equilibrium, 1e-6);
    }
    assert!(after[0].material.pressure_pa() > before[0].material.pressure_pa());
    assert!(after[1].material.pressure_pa() < before[1].material.pressure_pa());
    assert_eq!(audit.thermal.external_heat_j, 0.0);
    assert_eq!(audit.boundary_work_j, 0.0);
    near(audit.internal_energy_change_j, 0.0, 2e-6);
    near(
        audit.enthalpy_change_j,
        audit.pressure_enthalpy_change_j,
        2e-6,
    );
    assert!(audit.source_entropy_change_j_k > 0.0);
    assert!(audit.entropy.conduction_generation_j_k > 0.0);
    assert_eq!(audit.entropy.reservoir_entropy_input_j_k, 0.0);
    assert_eq!(next.thermal().mass_kg(), initial.thermal().mass_kg());
    assert_eq!(
        next.thermal().cells()[&CellId(0)].geometry,
        initial.thermal().cells()[&CellId(0)].geometry
    );
}

#[test]
fn bath_events_are_exact_and_refused_advances_preserve_sealed_state() {
    let law = law();
    let events = vec![
        birth(&law, 0, 650.0, 0.05, Some((750.0, 5.0))),
        TimedEvent {
            time_s: 10.125,
            event: Event::Boundary(CellId(0), boundary(Some((650.0, 5.0)))),
        },
    ];
    let initial = DenseWaterThermalNetwork::new(&law, events.clone(), budget(), || Ok(())).unwrap();
    let original = initial.profile(|| Ok(())).unwrap();
    let (whole, whole_audit) = initial
        .advance_to(25.125, controls(0.5), || Ok(()))
        .unwrap();
    let (warm, warm_audit) = initial
        .advance_to(10.125, controls(0.5), || Ok(()))
        .unwrap();
    assert_eq!(warm.time_s(), 10.125);
    assert_eq!(
        warm.thermal().cells()[&CellId(0)]
            .boundary
            .reservoir
            .unwrap()
            .temperature_k,
        650.0
    );
    let (cooled, cool_audit) = warm.advance_to(25.125, controls(0.5), || Ok(())).unwrap();
    assert_eq!(
        cooled.profile(|| Ok(())).unwrap(),
        whole.profile(|| Ok(())).unwrap()
    );
    assert!(warm_audit.thermal.external_heat_j > 0.0 && cool_audit.thermal.external_heat_j < 0.0);
    assert_eq!(whole_audit.thermal.events, 1);
    near(
        whole_audit.thermal.external_heat_j,
        warm_audit.thermal.external_heat_j + cool_audit.thermal.external_heat_j,
        1e-10,
    );
    near(
        whole_audit.entropy.reservoir_entropy_input_j_k,
        warm_audit.entropy.reservoir_entropy_input_j_k
            + cool_audit.entropy.reservoir_entropy_input_j_k,
        1e-10,
    );
    let mut calls = 0;
    let error = initial
        .advance_to(25.125, controls(0.5), || {
            calls += 1;
            if calls >= 200 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        })
        .err()
        .unwrap();
    assert_eq!(
        error,
        DenseWaterThermalError::Thermal(DepositionError::Cancelled)
    );
    let limited = DenseWaterThermalNetwork::new(
        &law,
        events,
        DepositionBudget {
            max_substeps: 1,
            ..budget()
        },
        || Ok(()),
    )
    .unwrap();
    assert!(matches!(
        limited.advance_to(1.0, controls(0.1), || Ok(())),
        Err(DenseWaterThermalError::Thermal(DepositionError::Budget(_)))
    ));
    let outside = DenseWaterThermalNetwork::new(
        &law,
        vec![birth(&law, 0, 650.0, 0.05, Some((1000.0, 5.0)))],
        budget(),
        || Ok(()),
    )
    .unwrap();
    assert!(matches!(
        outside.advance_to(100.0, controls(0.5), || Ok(())),
        Err(DenseWaterThermalError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(outside.time_s(), 0.0);
    assert_eq!(initial.profile(|| Ok(())).unwrap(), original);
    assert_eq!(limited.profile(|| Ok(())).unwrap(), original);
}
