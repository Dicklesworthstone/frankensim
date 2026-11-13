//! Analytic volume-loss oracles and the original thermal owner, without a
//! replacement acoustic or heat integrator. These are authored numerical
//! materials; they establish implementation behavior, not material calibration.

use super::*;
use fs_feec::acoustics::helmholtz::heterogeneous::{
    EquivalentFluidCoefficients, EquivalentFluidControls,
};
use fs_feec::acoustics::helmholtz::{HelmholtzFe, Medium};
use fs_rep_mesh::TetComplex;
use fs_thermal::enthalpy::{
    Boundary, DepositedCell, DepositionBudget, Reservoir, ThermalLink, TimedEvent,
};
use fs_thermal::observers::CoupledObservers;

mod walls;
mod feedback;
mod thermoviscous;

const IDS: [CellId; 2] = [CellId(91), CellId(7)];
const VOLUMES: [f64; 2] = [1.0 / 6.0, 1.0 / 3.0];
const CP: f64 = 200.0;
const FREQUENCY_HZ: f64 = 0.35;

#[derive(Clone, Copy)]
struct Caloric {
    maximum_h_j_kg: f64,
    identity: Option<[u8; 32]>,
    basis: ThermalStorageBasis,
}

impl ThermalLaw for Caloric {
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.basis
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        self.identity
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || !(-1000.0..=self.maximum_h_j_kg).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(300.0 + h / CP)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / CP
    }
}

fn caloric() -> Caloric {
    Caloric {
        maximum_h_j_kg: 1e6,
        identity: Some([31; 32]),
        basis: ThermalStorageBasis::DeclaredSensibleStorage,
    }
}

fn mesh() -> (TetComplex, Vec<[f64; 3]>) {
    (
        TetComplex::from_tets(5, vec![[0, 1, 2, 3], [0, 2, 1, 4]]),
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -2.0],
        ],
    )
}

fn geometry<'c>(mesh: &'c TetComplex, positions: &[[f64; 3]]) -> HelmholtzFe<'c> {
    HelmholtzFe::new(
        mesh,
        positions,
        1,
        Medium {
            density_kg_m3: 1.0,
            sound_speed_m_s: 2.0,
        },
    )
    .unwrap()
}

fn coefficients(element: usize) -> EquivalentFluidCoefficients {
    EquivalentFluidCoefficients {
        dynamic_density_kg_m3: [
            C64::new(1.0, -0.3),
            C64::new(2.0, -0.5),
        ][element],
        dynamic_bulk_modulus_pa: [
            C64::new(5.0, 0.8),
            C64::new(9.0, 1.4),
        ][element],
    }
}

fn prepared<'f, 'c>(geometry: &'f HelmholtzFe<'c>) -> PreparedEquivalentFluid<'f, 'c> {
    PreparedEquivalentFluid::new(
        geometry,
        FREQUENCY_HZ,
        &[],
        &[],
        EquivalentFluidControls::default(),
        &mut |element, _| Ok::<_, HeterogeneousError>(coefficients(element)),
        &mut || false,
    )
    .unwrap()
}

fn pressure(positions: &[[f64; 3]]) -> Vec<C64> {
    positions
        .iter()
        .map(|p| {
            C64::new(
                2.0 + 0.3 * p[0] - 0.2 * p[1] + 0.4 * p[2],
                -0.7 + 0.1 * p[0] + 0.6 * p[1] - 0.25 * p[2],
            )
        })
        .collect()
}

// Exact affine-field integral on a tetrahedron:
// integral |p|^2 = V (sum |p_i|^2 + |sum p_i|^2) / 20.
// No production quadrature, retained loss matrix, or acoustic power residual
// supplies this oracle. The prescribed field need not solve a boundary problem.
fn analytic_powers(mesh: &TetComplex, pressure: &[C64]) -> [[f64; 2]; 2] {
    let omega = core::f64::consts::TAU * FREQUENCY_HZ;
    let gradient_squared =
        0.3_f64.powi(2) + 0.2_f64.powi(2) + 0.4_f64.powi(2)
            + 0.1_f64.powi(2) + 0.6_f64.powi(2) + 0.25_f64.powi(2);
    core::array::from_fn(|element| {
        let values = mesh.tets[element].map(|id| pressure[id as usize]);
        let sum = values.iter().fold(C64::ZERO, |a, &b| a + b);
        let mean_square = (values.iter().map(|p| p.norm_sq()).sum::<f64>() + sum.norm_sq()) / 20.0;
        let material = coefficients(element);
        [
            material.dynamic_density_kg_m3.recip().im * VOLUMES[element]
                * gradient_squared / (2.0 * omega),
            -0.5 * omega * material.dynamic_bulk_modulus_pa.recip().im
                * VOLUMES[element] * mean_square,
        ]
    })
}

fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    }
}

fn schedule(sources_w: [f64; 2]) -> Vec<TimedEvent> {
    (0..2)
        .map(|element| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[element],
                segment_id: 100 + element as u64,
                body_id: 8,
                material_id: MaterialId(5),
                mass_kg: [0.5, 1.0][element],
                specific_h_j_kg: [20.0, 40.0][element],
                geometry: CellGeometry {
                    start_m: [1.0 / 3.0, 1.0 / 3.0, 0.0],
                    end_m: [0.0, 0.0, [1.0, -2.0][element]],
                    volume_m3: VOLUMES[element],
                    surface_area_m2: [(3.0 + 3.0_f64.sqrt()) / 2.0, 4.0][element],
                },
                boundary: boundary(sources_w[element]),
            }),
        })
        .collect()
}

fn config(max_dt_s: f64) -> StepConfig {
    StepConfig {
        max_dt_s,
        energy_tolerance_j: 1e-10,
    }
}

fn initial<'a>(schedule: Vec<TimedEvent>, law: &'a dyn ThermalLaw) -> DepositionSimulation<'a> {
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(5), law), (MaterialId(6), law)],
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 16,
            max_substeps: 2000,
        },
    )
    .unwrap()
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0
}

fn routes(law: &dyn ThermalLaw) -> [InteriorHeatCell<'_>; 2] {
    IDS.map(|cell| InteriorHeatCell { cell, caloric: law })
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance {tolerance:.3e}"
    );
}

#[test]
fn affine_element_losses_enter_original_cells_and_native_event_balance_once() {
    let (mesh, positions) = mesh();
    let geometry = geometry(&mesh, &positions);
    let prepared = prepared(&geometry);
    let pressure = pressure(&positions);
    let expected = analytic_powers(&mesh, &pressure);
    let law = caloric();
    let mut events = schedule([0.2, -0.05]);
    events.push(TimedEvent {
        time_s: 0.375,
        event: Event::Boundary(IDS[0], boundary(0.4)),
    });
    let thermal = initial(events, &law);
    let heating = InteriorAcousticHeating::try_new(
        &prepared, &pressure, &thermal, &routes(&law), 2, || Ok(()),
    )
    .unwrap();
    assert!(core::ptr::eq(heating.prepared(), &prepared));
    assert_eq!(heating.pressure_pa(), pressure.as_slice());
    assert_eq!(heating.cells().len(), 2);
    for (element, cell) in heating.cells().enumerate() {
        assert_eq!(cell.element, element);
        assert_eq!(cell.cell, IDS[element]);
        close(cell.volume_m3, VOLUMES[element], 1e-14);
        close(cell.density_dissipation_w, expected[element][0], 2e-13);
        close(cell.compressional_dissipation_w, expected[element][1], 2e-13);
        assert!(cell.density_dissipation_w > 0.0 && cell.compressional_dissipation_w > 0.0);
    }
    let step = heating.advance_to(&thermal, 1.0, config(0.13), || Ok(())).unwrap();
    let authored_j = [0.2 * 0.375 + 0.4 * 0.625, -0.05];
    for element in 0..2 {
        let acoustic_j = expected[element].iter().sum::<f64>();
        let original = &thermal.cells()[&IDS[element]];
        let final_cell = &step.thermal.cells()[&IDS[element]];
        assert_eq!(step.cell_energy_j[element].0, IDS[element]);
        close(step.cell_energy_j[element].1, acoustic_j, 3e-13);
        close(
            (final_cell.specific_h_j_kg - original.specific_h_j_kg) * original.mass_kg,
            acoustic_j + authored_j[element],
            2e-12,
        );
        close(
            step.thermal.temperature_k(IDS[element]).unwrap(),
            300.0 + (original.specific_h_j_kg
                + (acoustic_j + authored_j[element]) / original.mass_kg) / CP,
            2e-12,
        );
        assert_eq!(final_cell.geometry, original.geometry);
        assert_eq!(final_cell.mass_kg, original.mass_kg);
    }
    let density_j = expected.iter().map(|p| p[0]).sum::<f64>();
    let compression_j = expected.iter().map(|p| p[1]).sum::<f64>();
    close(step.acoustic.density_dissipation_j, density_j, 3e-13);
    close(step.acoustic.compressional_dissipation_j, compression_j, 3e-13);
    close(step.acoustic.deposited_j, density_j + compression_j, 3e-13);
    close(
        step.heat_report.external_heat_j,
        step.acoustic.deposited_j + authored_j.iter().sum::<f64>(),
        2e-12,
    );
    close(step.heat_report.residual_j, 0.0, 2e-12);
    assert_eq!(step.heat_report.events, 1);
    assert!(step.heat_report.substeps > 1);
    assert_eq!(step.thermal.cells()[&IDS[0]].boundary.source_w, 0.4);
    assert_eq!(step.thermal.cells()[&IDS[1]].boundary.source_w, -0.05);
    assert_eq!(thermal.cells()[&IDS[0]].boundary.source_w, 0.2);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn composed_observer_shares_native_cfl_contacts_reservoirs_and_split_history() {
    let (mesh, positions) = mesh();
    let geometry = geometry(&mesh, &positions);
    let prepared = prepared(&geometry);
    let pressure = pressure(&positions);
    let expected = analytic_powers(&mesh, &pressure);
    let law = caloric();
    let mut events = schedule([0.2, -0.05]);
    if let Event::Birth(cell) = &mut events[1].event {
        cell.boundary.reservoir = Some(Reservoir {
            area_m2: 0.25,
            conductance_w_k: 320.0,
            temperature_k: 299.9,
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 3, a: IDS[0], b: IDS[1], area_m2: 0.5, conductance_w_k: 640.0,
        }),
    });
    events.push(TimedEvent {
        time_s: 0.375,
        event: Event::Boundary(IDS[0], boundary(0.4)),
    });
    let thermal = initial(events.clone(), &law);
    let source = InteriorAcousticHeating::try_new(
        &prepared, &pressure, &thermal, &routes(&law), 2, || Ok(()),
    )
    .unwrap();
    let observers = CoupledObservers::new((), source.observer());
    // The native CFL bound is exactly 0.125 s. The event at 0.375 s and
    // the split at 0.5 s preserve the same accepted step grid as the whole run.
    let (middle, middle_observers, first_report) =
        thermal.advance_observed(&observers, 0.5, config(10.0), || Ok(())).unwrap();
    let (actual, accepted_observers, second_report) = middle
        .advance_observed(&middle_observers, 1.0, config(10.0), || Ok(()))
        .unwrap();
    let (whole, whole_observers, _) =
        thermal.advance_observed(&observers, 1.0, config(10.0), || Ok(())).unwrap();

    // Independent reference: author the analytically integrated acoustic
    // power directly on the SAME native thermal owner, including its events.
    let extra = expected.map(|p| p[0] + p[1]);
    for event in &mut events {
        match &mut event.event {
            Event::Birth(cell) => {
                let element = IDS.iter().position(|&id| id == cell.id).unwrap();
                cell.boundary.source_w += extra[element];
            }
            Event::Boundary(id, boundary) => {
                let element = IDS.iter().position(|&candidate| candidate == *id).unwrap();
                boundary.source_w += extra[element];
            }
            _ => {}
        }
    }
    let reference = initial(events, &law);
    let (reference_middle, native_first) =
        reference.advance_to(0.5, config(10.0), || Ok(())).unwrap();
    let (reference_end, native_second) =
        reference_middle.advance_to(1.0, config(10.0), || Ok(())).unwrap();
    assert!(first_report.substeps + second_report.substeps > 2);
    assert_eq!(first_report.substeps, native_first.substeps);
    assert_eq!(second_report.substeps, native_second.substeps);
    for id in IDS {
        close(
            actual.cells()[&id].specific_h_j_kg,
            reference_end.cells()[&id].specific_h_j_kg,
            2e-10,
        );
        close(
            actual.cells()[&id].specific_h_j_kg,
            whole.cells()[&id].specific_h_j_kg,
            2e-10,
        );
        assert_eq!(actual.cells()[&id].geometry, thermal.cells()[&id].geometry);
    }
    assert_eq!(actual.links(), thermal.links());
    assert_eq!(actual.cells()[&IDS[0]].boundary.source_w, 0.4);
    assert_eq!(actual.cells()[&IDS[1]].boundary, thermal.cells()[&IDS[1]].boundary);
    close(
        first_report.external_heat_j + second_report.external_heat_j,
        native_first.external_heat_j + native_second.external_heat_j,
        2e-10,
    );
    let expected_j = expected.iter().flatten().sum::<f64>();
    close(accepted_observers.second.energy().unwrap().deposited_j, expected_j, 3e-13);
    close(whole_observers.second.energy().unwrap().deposited_j, expected_j, 3e-13);
    close(middle_observers.second.energy().unwrap().deposited_j, 0.5 * expected_j, 3e-13);
    assert_eq!(observers.second.energy().unwrap().deposited_j, 0.0);

    let zero_pressure = vec![C64::ZERO; geometry.ndof()];
    let zero = InteriorAcousticHeating::try_new(
        &prepared, &zero_pressure, &thermal, &routes(&law), 2, || Ok(()),
    )
    .unwrap();
    let zero_step = zero.advance_to(&thermal, 1.0, config(10.0), || Ok(())).unwrap();
    let native = thermal.advance_to(1.0, config(10.0), || Ok(())).unwrap();
    assert_eq!(zero_step.thermal.cells(), native.0.cells());
    assert_eq!(zero_step.heat_report, native.1);
    assert_eq!(zero_step.acoustic.deposited_j, 0.0);
}

#[test]
fn original_volume_caloric_object_and_cell_binding_cannot_be_substituted() {
    let (mesh, positions) = mesh();
    let geometry = geometry(&mesh, &positions);
    let prepared = prepared(&geometry);
    let pressure = pressure(&positions);
    let law = caloric();
    let twin = law;
    let events = schedule([0.0; 2]);
    let thermal = initial(events.clone(), &law);
    let mapping = routes(&law);
    let heating = InteriorAcousticHeating::try_new(
        &prepared, &pressure, &thermal, &mapping, 2, || Ok(()),
    )
    .unwrap();
    for cap in [0, 1, MAX_INTERIOR_HEATING_CELLS + 1] {
        assert!(matches!(
            InteriorAcousticHeating::try_new(
                &prepared, &pressure, &thermal, &mapping, cap, || Ok(()),
            ),
            Err(InteriorHeatingError::Thermal(DepositionError::Budget(_)))
        ));
    }
    assert!(InteriorAcousticHeating::try_new(
        &prepared, &pressure, &thermal, &mapping[..1], 2, || Ok(()),
    ).is_err());
    for bad in [
        [mapping[0], mapping[0]],
        [mapping[1], mapping[0]],
        [mapping[0], InteriorHeatCell { cell: CellId(999), caloric: &law }],
        routes(&twin),
    ] {
        assert!(matches!(
            InteriorAcousticHeating::try_new(
                &prepared, &pressure, &thermal, &bad, 2, || Ok(()),
            ),
            Err(InteriorHeatingError::Invalid(_))
        ));
    }
    // Equal declared identity AND equal temperature do not make a new object
    // the original caloric law borrowed by this timeline.
    assert_eq!(law.caloric_identity(), twin.caloric_identity());
    assert_eq!(law.temperature_k(20.0), twin.temperature_k(20.0));
    let substituted = initial(events.clone(), &twin);
    assert!(matches!(
        heating.advance_to(&substituted, 0.1, config(0.1), || Ok(())),
        Err(InteriorHeatingError::Invalid(_))
    ));
    for mutation in 0..6 {
        let mut changed = events.clone();
        let Event::Birth(cell) = &mut changed[0].event else { panic!("birth") };
        match mutation {
            0 => cell.geometry.volume_m3 *= 1.01,
            1 => cell.geometry.start_m[0] += 0.01,
            2 => cell.mass_kg *= 1.01,
            3 => cell.body_id += 1,
            4 => cell.segment_id += 1,
            _ => cell.material_id = MaterialId(6),
        }
        let changed = initial(changed, &law);
        assert!(matches!(
            heating.advance_to(&changed, 0.1, config(0.1), || Ok(())),
            Err(InteriorHeatingError::Invalid(_))
        ));
        if mutation == 0 {
            assert!(matches!(
                InteriorAcousticHeating::try_new(
                    &prepared, &pressure, &changed, &mapping, 2, || Ok(()),
                ),
                Err(InteriorHeatingError::Invalid(_))
            ));
        }
    }
    assert!(matches!(
        InteriorAcousticHeating::try_new(
            &prepared, &pressure[..pressure.len() - 1], &thermal, &mapping, 2, || Ok(()),
        ),
        Err(InteriorHeatingError::Acoustic(_))
    ));
    let mut nonfinite = pressure.clone();
    nonfinite[1] = C64::new(f64::NAN, 0.0);
    assert!(matches!(
        InteriorAcousticHeating::try_new(
            &prepared, &nonfinite, &thermal, &mapping, 2, || Ok(()),
        ),
        Err(InteriorHeatingError::Acoustic(_))
    ));

    // Heat-only coupling preserves an absent declared caloric identity while
    // still requiring the exact original object and its storage basis.
    let anonymous = Caloric { identity: None, basis: ThermalStorageBasis::InternalEnergy, ..law };
    let anonymous_thermal = initial(events, &anonymous);
    let anonymous_heat = InteriorAcousticHeating::try_new(
        &prepared, &pressure, &anonymous_thermal, &routes(&anonymous), 2, || Ok(()),
    )
    .unwrap();
    let anonymous_step = anonymous_heat
        .advance_to(&anonymous_thermal, 0.1, config(0.1), || Ok(()))
        .unwrap();
    assert_eq!(
        anonymous_step.thermal.material_storage_basis(MaterialId(5)).unwrap(),
        ThermalStorageBasis::InternalEnergy,
    );
    assert_eq!(anonymous_step.thermal.material_caloric_identity(MaterialId(5)).unwrap(), None);
}

#[test]
fn cancellation_material_domain_and_endpoint_events_roll_back_both_owners() {
    let (mesh, positions) = mesh();
    let geometry = geometry(&mesh, &positions);
    let prepared = prepared(&geometry);
    let pressure = pressure(&positions);
    let law = caloric();
    let events = schedule([0.2, 0.1]);
    let thermal = initial(events.clone(), &law);
    let original_cells = thermal.cells().clone();
    let mut preparation_polls = 0;
    let heating = InteriorAcousticHeating::try_new(
        &prepared, &pressure, &thermal, &routes(&law), 2, || {
            preparation_polls += 1;
            Ok(())
        },
    )
    .unwrap();
    let mut polls = 0;
    assert!(matches!(
        InteriorAcousticHeating::try_new(
            &prepared, &pressure, &thermal, &routes(&law), 2, || {
                polls += 1;
                if polls >= preparation_polls / 2 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        ),
        Err(InteriorHeatingError::Thermal(DepositionError::Cancelled))
    ));
    let observer = heating.observer();
    let mut advance_polls = 0;
    let (_, accepted, _) = thermal
        .advance_observed(&observer, 1.0, config(0.1), || {
            advance_polls += 1;
            Ok(())
        })
        .unwrap();
    assert!(accepted.energy().unwrap().deposited_j > 0.0);
    let mut polls = 0;
    assert!(matches!(
        thermal.advance_observed(&observer, 1.0, config(0.1), || {
            polls += 1;
            if polls >= advance_polls * 3 / 4 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(InteriorHeatingError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(observer.energy().unwrap().deposited_j, 0.0);
    assert!(observer.cell_energy_j().unwrap().iter().all(|&(_, energy)| energy == 0.0));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells(), &original_cells);

    // This event fails only after native substeps have reached the requested
    // endpoint. The original event cursor, cells and observer account survive.
    let mut failing_events = events;
    failing_events.push(TimedEvent {
        time_s: 0.5,
        event: Event::Boundary(CellId(999), boundary(0.0)),
    });
    let event_thermal = initial(failing_events, &law);
    let event_heat = InteriorAcousticHeating::try_new(
        &prepared, &pressure, &event_thermal, &routes(&law), 2, || Ok(()),
    )
    .unwrap();
    let event_observer = event_heat.observer();
    assert!(matches!(
        event_thermal.advance_observed(&event_observer, 0.5, config(0.1), || Ok(())),
        Err(InteriorHeatingError::Thermal(DepositionError::InvalidInput(_)))
    ));
    assert_eq!(event_thermal.time_s(), 0.0);
    assert_eq!(event_thermal.cells(), &original_cells);
    assert_eq!(event_observer.energy().unwrap().deposited_j, 0.0);

    let limited = Caloric { maximum_h_j_kg: 0.0001, ..law };
    let mut limited_events = schedule([0.2, 0.1]);
    for event in &mut limited_events {
        if let Event::Birth(cell) = &mut event.event {
            cell.specific_h_j_kg = 0.0;
        }
    }
    let limited_thermal = initial(limited_events, &limited);
    let limited_cells = limited_thermal.cells().clone();
    let limited_heat = InteriorAcousticHeating::try_new(
        &prepared, &pressure, &limited_thermal, &routes(&limited), 2, || Ok(()),
    )
    .unwrap();
    let limited_observer = limited_heat.observer();
    assert!(matches!(
        limited_thermal.advance_observed(&limited_observer, 0.5, config(0.1), || Ok(())),
        Err(InteriorHeatingError::Thermal(DepositionError::MaterialDomain))
    ));
    assert_eq!(limited_thermal.cells(), &limited_cells);
    assert_eq!(limited_thermal.time_s(), 0.0);
    assert_eq!(limited_observer.energy().unwrap().deposited_j, 0.0);
}
