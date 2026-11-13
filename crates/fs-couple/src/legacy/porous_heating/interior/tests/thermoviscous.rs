//! Independent polynomial loss moments and the original thermal timeline.
//!
//! These manufactured full P2 fields test loss observation and heat transfer;
//! they are not claimed to solve the mixed acoustic boundary-value problem.

use super::*;
use crate::porous_heating::interior::thermoviscous::{
    ThermoviscousHeatingEnergy, ThermoviscousHeatingError, ThermoviscousViscousHeating,
};
use fs_feec::acoustics::helmholtz::multiport::BoundaryMode;
use fs_feec::acoustics::thermoviscous::{
    ThermoviscousBoundary, ThermoviscousElementPower, ThermoviscousError, ThermoviscousFe,
    ThermoviscousLimits, ThermoviscousMedium,
};
use fs_thermal::observers::ObserverCouplingError;

const LNS_FREQUENCY_HZ: f64 = 137.0;
const COMPLEX_CP: C64 = C64::new(970.0, -43.0);
const VELOCITY_AMPLITUDE: C64 = C64::new(0.7, -0.4);
const TEMPERATURE_AMPLITUDE: C64 = C64::new(-0.2, 0.6);

fn medium() -> ThermoviscousMedium {
    ThermoviscousMedium {
        // The original unequal-volume thermal fixture has masses 0.5 and 1 kg.
        density_kg_m3: 3.0,
        pressure_pa: 3.0 * 287.05 * 300.0,
        temperature_k: 300.0,
        cp_j_kg_k: 1005.0,
        dynamic_viscosity_pa_s: 1.81e-5,
        thermal_conductivity_w_m_k: 0.0257,
        bulk_viscosity_pa_s: 2.7e-5,
    }
}

fn lns<'c>(mesh: &'c TetComplex, positions: &[[f64; 3]]) -> ThermoviscousFe<'c> {
    let boundary = ThermoviscousBoundary {
        no_slip_faces: vec![],
        isothermal_faces: vec![],
        ports: vec![BoundaryMode::uniform(
            &geometry(mesh, positions).boundary_faces(),
        )],
    };
    ThermoviscousFe::new(
        mesh,
        positions,
        medium(),
        &boundary,
        ThermoviscousLimits {
            max_dofs: 128,
            max_dense_bytes: 48 * 128 * 128,
            max_elements: 2,
            max_element_aspect_ratio: 100.0,
        },
        &mut || false,
    )
    .unwrap()
}

fn affine_fields(mesh: &TetComplex, positions: &[[f64; 3]]) -> (Vec<[C64; 3]>, Vec<C64>) {
    // Native P2 is hierarchical: an affine field has zero edge coefficients.
    let count = positions.len() + mesh.edges.len();
    let mut velocity = vec![[C64::ZERO; 3]; count];
    let mut temperature = vec![C64::ZERO; count];
    for (i, &[x, y, z]) in positions.iter().enumerate() {
        velocity[i] = [
            VELOCITY_AMPLITUDE.scale(2.0 * x - y),
            VELOCITY_AMPLITUDE.scale(3.0 * y + z),
            VELOCITY_AMPLITUDE.scale(-x + 4.0 * z),
        ];
        temperature[i] = TEMPERATURE_AMPLITUDE.scale(1.0 + 2.0 * x - 3.0 * y + z);
    }
    (velocity, temperature)
}

fn affine_powers() -> [ThermoviscousElementPower; 2] {
    let medium = medium();
    std::array::from_fn(|element| {
        let volume = VOLUMES[element];
        // div(v)=9a; dev(sym(grad v)) has squared norm 7|a|²/2.
        // |grad T'|²=14|b|². Exact simplex moments of
        // (1+2x-3y+z)² are 17/10 and 4/5 on these original tets.
        let thermal = 7.0 * medium.thermal_conductivity_w_m_k / medium.temperature_k
            * volume
            * TEMPERATURE_AMPLITUDE.norm_sq();
        ThermoviscousElementPower {
            volume_m3: volume,
            mass_kg: medium.density_kg_m3 * volume,
            shear_heating_w: 3.5
                * medium.dynamic_viscosity_pa_s
                * volume
                * VELOCITY_AMPLITUDE.norm_sq(),
            bulk_heating_w: 40.5
                * medium.bulk_viscosity_pa_s
                * volume
                * VELOCITY_AMPLITUDE.norm_sq(),
            conductive_entropy_w_k: thermal / medium.temperature_k,
            thermal_availability_loss_w: thermal,
            caloric_availability_loss_w: -0.5
                * std::f64::consts::TAU
                * LNS_FREQUENCY_HZ
                * medium.density_kg_m3
                * COMPLEX_CP.im
                / medium.temperature_k
                * volume
                * [1.7, 0.8][element]
                * TEMPERATURE_AMPLITUDE.norm_sq(),
        }
    })
}

fn relative(actual: f64, expected: f64) {
    close(actual, expected, 4e-12 * expected.abs().max(1e-20));
}

fn check_energy(
    energy: ThermoviscousHeatingEnergy,
    powers: &[ThermoviscousElementPower],
    duration_s: f64,
) {
    let expected = |value: fn(&ThermoviscousElementPower) -> f64| {
        duration_s * powers.iter().map(value).sum::<f64>()
    };
    relative(energy.shear_heating_j, expected(|p| p.shear_heating_w));
    relative(energy.bulk_heating_j, expected(|p| p.bulk_heating_w));
    relative(
        energy.conductive_entropy_j_k,
        expected(|p| p.conductive_entropy_w_k),
    );
    relative(
        energy.thermal_availability_loss_j,
        expected(|p| p.thermal_availability_loss_w),
    );
    relative(
        energy.caloric_availability_loss_j,
        expected(|p| p.caloric_availability_loss_w),
    );
    relative(
        energy.deposited_j,
        expected(|p| p.shear_heating_w + p.bulk_heating_w),
    );
}

#[test]
fn affine_lns_positive_forms_heat_native_cells_once_across_events() {
    let (mesh, positions) = mesh();
    let model = lns(&mesh, &positions);
    let prepared = model
        .prepare_with_heat_capacity(LNS_FREQUENCY_HZ, COMPLEX_CP, &mut || false)
        .unwrap();
    let (velocity, temperature) = affine_fields(&mesh, &positions);
    let expected = affine_powers();
    let law = caloric();
    let mut events = schedule([0.2, -0.05]);
    events.push(TimedEvent {
        time_s: 0.375,
        event: Event::Boundary(IDS[0], boundary(0.4)),
    });
    let thermal = initial(events, &law);
    let heating = ThermoviscousViscousHeating::try_new(
        &prepared,
        &velocity,
        &temperature,
        &thermal,
        &routes(&law),
        2,
        || Ok(()),
    )
    .unwrap();
    assert!(core::ptr::eq(heating.prepared(), &prepared));
    assert_eq!(heating.velocity_m_s(), velocity.as_slice());
    assert_eq!(heating.temperature_k(), temperature.as_slice());
    assert_eq!(heating.cells().len(), 2);
    for (element, cell) in heating.cells().enumerate() {
        assert_eq!((cell.element, cell.cell), (element, IDS[element]));
        let power = cell.power;
        let expected = expected[element];
        relative(power.volume_m3, expected.volume_m3);
        relative(power.mass_kg, expected.mass_kg);
        relative(power.shear_heating_w, expected.shear_heating_w);
        relative(power.bulk_heating_w, expected.bulk_heating_w);
        relative(
            power.conductive_entropy_w_k,
            expected.conductive_entropy_w_k,
        );
        relative(
            power.thermal_availability_loss_w,
            expected.thermal_availability_loss_w,
        );
        relative(
            power.caloric_availability_loss_w,
            expected.caloric_availability_loss_w,
        );
    }
    let step = heating
        .advance_to(&thermal, 1.0, config(0.13), || Ok(()))
        .unwrap();
    let authored_j = [0.2 * 0.375 + 0.4 * 0.625, -0.05];
    for element in 0..2 {
        let heat_j = expected[element].viscous_heat_w();
        assert_eq!(step.cell_energy_j[element].0, IDS[element]);
        relative(step.cell_energy_j[element].1, heat_j);
        let old = &thermal.cells()[&IDS[element]];
        let new = &step.thermal.cells()[&IDS[element]];
        close(
            old.mass_kg * (new.specific_h_j_kg - old.specific_h_j_kg),
            authored_j[element] + heat_j,
            2e-12,
        );
        assert_eq!(new.geometry, old.geometry);
        assert_eq!(new.mass_kg, old.mass_kg);
    }
    check_energy(step.acoustic, &expected, 1.0);
    assert!(step.acoustic.caloric_availability_loss_j > step.acoustic.deposited_j);
    close(
        step.heat_report.external_heat_j,
        authored_j.iter().sum::<f64>() + step.acoustic.deposited_j,
        2e-12,
    );
    close(step.heat_report.residual_j, 0.0, 2e-12);
    assert_eq!(step.heat_report.events, 1);
    assert!(step.heat_report.substeps > 1);
    assert_eq!(step.thermal.cells()[&IDS[0]].boundary.source_w, 0.4);
    assert_eq!(thermal.cells()[&IDS[0]].boundary.source_w, 0.2);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn quadratic_temperature_availability_is_observed_without_caloric_deposition() {
    let (mesh, positions) = mesh();
    let model = lns(&mesh, &positions);
    let prepared = model
        .prepare_with_heat_capacity(LNS_FREQUENCY_HZ, COMPLEX_CP, &mut || false)
        .unwrap();
    // Native edge bubbles are lambda_a*lambda_b, so x² has the edge
    // coefficient -(x_a-x_b)², not its quadratic midpoint nodal value.
    let mut temperature = positions
        .iter()
        .map(|p| TEMPERATURE_AMPLITUDE.scale(p[0] * p[0]))
        .collect::<Vec<_>>();
    for &[a, b] in &mesh.edges {
        let dx = positions[a as usize][0] - positions[b as usize][0];
        temperature.push(TEMPERATURE_AMPLITUDE.scale(-dx * dx));
    }
    let velocity = vec![[C64::ZERO; 3]; temperature.len()];
    let law = caloric();
    let thermal = initial(schedule([0.2, -0.05]), &law);
    let heating = ThermoviscousViscousHeating::try_new(
        &prepared,
        &velocity,
        &temperature,
        &thermal,
        &routes(&law),
        2,
        || Ok(()),
    )
    .unwrap();
    let medium = medium();
    let expected = VOLUMES.map(|volume| {
        // On both physical tetrahedra, mean x²=1/10 and mean x⁴=1/35.
        let thermal = medium.thermal_conductivity_w_m_k * volume / (5.0 * medium.temperature_k)
            * TEMPERATURE_AMPLITUDE.norm_sq();
        ThermoviscousElementPower {
            volume_m3: volume,
            mass_kg: medium.density_kg_m3 * volume,
            shear_heating_w: 0.0,
            bulk_heating_w: 0.0,
            conductive_entropy_w_k: thermal / medium.temperature_k,
            thermal_availability_loss_w: thermal,
            caloric_availability_loss_w: -0.5
                * std::f64::consts::TAU
                * LNS_FREQUENCY_HZ
                * medium.density_kg_m3
                * COMPLEX_CP.im
                / medium.temperature_k
                * volume
                / 35.0
                * TEMPERATURE_AMPLITUDE.norm_sq(),
        }
    });
    let step = heating
        .advance_to(&thermal, 1.0, config(0.13), || Ok(()))
        .unwrap();
    let native = thermal.advance_to(1.0, config(0.13), || Ok(())).unwrap();
    assert_eq!(step.thermal.cells(), native.0.cells());
    assert_eq!(step.heat_report, native.1);
    assert_eq!(step.acoustic.deposited_j, 0.0);
    assert!(step.cell_energy_j.iter().all(|&(_, heat)| heat == 0.0));
    assert!(step.acoustic.thermal_availability_loss_j > 0.0);
    assert!(step.acoustic.caloric_availability_loss_j > 0.0);
    check_energy(step.acoustic, &expected, 1.0);
}

#[test]
fn composed_continuation_uses_native_cfl_and_cancellation_preserves_both_histories() {
    let (mesh, positions) = mesh();
    let model = lns(&mesh, &positions);
    let prepared = model
        .prepare_with_heat_capacity(LNS_FREQUENCY_HZ, COMPLEX_CP, &mut || false)
        .unwrap();
    let (velocity, temperature) = affine_fields(&mesh, &positions);
    let expected = affine_powers();
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
            id: 3,
            a: IDS[0],
            b: IDS[1],
            area_m2: 0.5,
            conductance_w_k: 640.0,
        }),
    });
    events.push(TimedEvent {
        time_s: 0.375,
        event: Event::Boundary(IDS[0], boundary(0.4)),
    });
    let thermal = initial(events.clone(), &law);
    let heating = ThermoviscousViscousHeating::try_new(
        &prepared,
        &velocity,
        &temperature,
        &thermal,
        &routes(&law),
        2,
        || Ok(()),
    )
    .unwrap();
    let observers = CoupledObservers::new((), heating.observer());
    // Original masses, Cp and conductances give a 0.125 s native CFL bound.
    // Both the event and split occur on the same grid as the whole advance.
    let (middle, middle_observers, first_report) = thermal
        .advance_observed(&observers, 0.5, config(10.0), || Ok(()))
        .unwrap();
    let (actual, accepted, second_report) = middle
        .advance_observed(&middle_observers, 1.0, config(10.0), || Ok(()))
        .unwrap();
    let mut whole_polls = 0;
    let (whole, whole_observers, _) = thermal
        .advance_observed(&observers, 1.0, config(10.0), || {
            whole_polls += 1;
            Ok(())
        })
        .unwrap();

    // Independent source injection on the SAME native transport owner.
    // The polynomial oracle supplies only shear+bulk, including after events.
    for event in &mut events {
        match &mut event.event {
            Event::Birth(cell) => {
                let element = IDS.iter().position(|&id| id == cell.id).unwrap();
                cell.boundary.source_w += expected[element].viscous_heat_w();
            }
            Event::Boundary(id, boundary) => {
                let element = IDS.iter().position(|&candidate| candidate == *id).unwrap();
                boundary.source_w += expected[element].viscous_heat_w();
            }
            _ => {}
        }
    }
    let reference = initial(events, &law);
    let (reference_middle, native_first) =
        reference.advance_to(0.5, config(10.0), || Ok(())).unwrap();
    let (reference_end, native_second) = reference_middle
        .advance_to(1.0, config(10.0), || Ok(()))
        .unwrap();
    assert_eq!(first_report.substeps, 4);
    assert_eq!(second_report.substeps, 4);
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
    check_energy(accepted.second.energy().unwrap(), &expected, 1.0);
    check_energy(whole_observers.second.energy().unwrap(), &expected, 1.0);
    check_energy(middle_observers.second.energy().unwrap(), &expected, 0.5);
    for (element, &(id, heat)) in accepted.second.cell_energy_j().unwrap().iter().enumerate() {
        assert_eq!(id, IDS[element]);
        relative(heat, expected[element].viscous_heat_w());
    }
    close(
        first_report.external_heat_j + second_report.external_heat_j,
        native_first.external_heat_j + native_second.external_heat_j,
        2e-10,
    );

    let original_cells = thermal.cells().clone();
    let mut polls = 0;
    assert!(whole_polls > 16);
    assert!(matches!(
        thermal.advance_observed(&observers, 1.0, config(10.0), || {
            polls += 1;
            if polls >= whole_polls * 3 / 4 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(ObserverCouplingError::Thermal(DepositionError::Cancelled))
            | Err(ObserverCouplingError::Second(
                ThermoviscousHeatingError::Thermal(DepositionError::Cancelled)
            ))
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells(), &original_cells);
    check_energy(observers.second.energy().unwrap(), &expected, 0.0);
    assert!(
        observers
            .second
            .cell_energy_j()
            .unwrap()
            .iter()
            .all(|&(_, heat)| heat == 0.0)
    );
    let (retry, retried_observers, _) = thermal
        .advance_observed(&observers, 1.0, config(10.0), || Ok(()))
        .unwrap();
    assert_eq!(retry.cells(), whole.cells());
    assert_eq!(
        retried_observers.second.energy().unwrap(),
        whole_observers.second.energy().unwrap(),
    );
}

#[test]
fn original_mass_volume_fields_and_caloric_binding_are_required_before_publication() {
    let (mesh, positions) = mesh();
    let model = lns(&mesh, &positions);
    let prepared = model
        .prepare_with_heat_capacity(LNS_FREQUENCY_HZ, COMPLEX_CP, &mut || false)
        .unwrap();
    let (velocity, temperature) = affine_fields(&mesh, &positions);
    let law = caloric();
    let twin = law;
    let events = schedule([0.0; 2]);
    let thermal = initial(events.clone(), &law);
    let mapping = routes(&law);
    let mut preparation_polls = 0;
    let heating = ThermoviscousViscousHeating::try_new(
        &prepared,
        &velocity,
        &temperature,
        &thermal,
        &mapping,
        2,
        || {
            preparation_polls += 1;
            Ok(())
        },
    )
    .unwrap();
    // Cell admission precedes even malformed-field quadrature.
    for cap in [0, 1, MAX_INTERIOR_HEATING_CELLS + 1] {
        assert!(matches!(
            ThermoviscousViscousHeating::try_new(
                &prepared,
                &[],
                &[],
                &thermal,
                &mapping,
                cap,
                || Ok(()),
            ),
            Err(ThermoviscousHeatingError::Thermal(DepositionError::Budget(
                _
            )))
        ));
    }
    for bad in [
        [mapping[0], mapping[0]],
        [mapping[1], mapping[0]],
        [
            mapping[0],
            InteriorHeatCell {
                cell: CellId(999),
                caloric: &law,
            },
        ],
        routes(&twin),
    ] {
        assert!(matches!(
            ThermoviscousViscousHeating::try_new(
                &prepared,
                &velocity,
                &temperature,
                &thermal,
                &bad,
                2,
                || Ok(()),
            ),
            Err(ThermoviscousHeatingError::Binding(
                InteriorHeatingError::Invalid(_)
            ))
        ));
    }
    assert!(matches!(
        ThermoviscousViscousHeating::try_new(
            &prepared,
            &velocity,
            &temperature,
            &thermal,
            &mapping[..1],
            2,
            || Ok(()),
        ),
        Err(ThermoviscousHeatingError::Binding(
            InteriorHeatingError::Invalid(_)
        ))
    ));
    for mutation in 0..6 {
        let mut changed = events.clone();
        let Event::Birth(cell) = &mut changed[0].event else {
            panic!("birth")
        };
        match mutation {
            0 => cell.geometry.volume_m3 *= 1.01,
            1 => cell.mass_kg *= 1.01,
            2 => cell.geometry.start_m[0] += 0.01,
            3 => cell.body_id += 1,
            4 => cell.segment_id += 1,
            _ => cell.material_id = MaterialId(6),
        }
        let changed = initial(changed, &law);
        if mutation < 2 {
            assert!(matches!(
                ThermoviscousViscousHeating::try_new(
                    &prepared,
                    &velocity,
                    &temperature,
                    &changed,
                    &mapping,
                    2,
                    || Ok(()),
                ),
                Err(ThermoviscousHeatingError::Binding(
                    InteriorHeatingError::Invalid(_)
                ))
            ));
        }
        assert!(matches!(
            heating.advance_to(&changed, 0.1, config(0.1), || Ok(())),
            Err(ThermoviscousHeatingError::Binding(
                InteriorHeatingError::Invalid(_)
            ))
        ));
    }
    let substituted = initial(events, &twin);
    assert!(matches!(
        heating.advance_to(&substituted, 0.1, config(0.1), || Ok(())),
        Err(ThermoviscousHeatingError::Binding(
            InteriorHeatingError::Invalid(_)
        ))
    ));
    for (v, t) in [
        (&velocity[..velocity.len() - 1], temperature.as_slice()),
        (velocity.as_slice(), &temperature[..temperature.len() - 1]),
    ] {
        assert!(matches!(
            ThermoviscousViscousHeating::try_new(&prepared, v, t, &thermal, &mapping, 2, || Ok(()),),
            Err(ThermoviscousHeatingError::Acoustic(
                ThermoviscousError::InvalidInput(_)
            ))
        ));
    }
    let mut bad_velocity = velocity.clone();
    bad_velocity[0][1] = C64::new(f64::INFINITY, 0.0);
    let mut bad_temperature = temperature.clone();
    bad_temperature[1] = C64::new(0.0, f64::NAN);
    for (v, t) in [
        (bad_velocity.as_slice(), temperature.as_slice()),
        (velocity.as_slice(), bad_temperature.as_slice()),
    ] {
        assert!(matches!(
            ThermoviscousViscousHeating::try_new(&prepared, v, t, &thermal, &mapping, 2, || Ok(()),),
            Err(ThermoviscousHeatingError::Acoustic(
                ThermoviscousError::NonFinite
            ))
        ));
    }
    let mut polls = 0;
    assert!(matches!(
        ThermoviscousViscousHeating::try_new(
            &prepared,
            &velocity,
            &temperature,
            &thermal,
            &mapping,
            2,
            || {
                polls += 1;
                if polls >= preparation_polls / 2 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        ),
        Err(ThermoviscousHeatingError::Thermal(
            DepositionError::Cancelled
        ))
    ));

    let limited = Caloric {
        maximum_h_j_kg: 1e-8,
        ..law
    };
    let mut limited_events = schedule([0.0; 2]);
    for event in &mut limited_events {
        if let Event::Birth(cell) = &mut event.event {
            cell.specific_h_j_kg = 0.0;
        }
    }
    let limited_thermal = initial(limited_events, &limited);
    let limited_cells = limited_thermal.cells().clone();
    let limited_heat = ThermoviscousViscousHeating::try_new(
        &prepared,
        &velocity,
        &temperature,
        &limited_thermal,
        &routes(&limited),
        2,
        || Ok(()),
    )
    .unwrap();
    let observer = limited_heat.observer();
    assert!(matches!(
        limited_thermal.advance_observed(&observer, 0.5, config(0.1), || Ok(())),
        Err(ThermoviscousHeatingError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
    assert_eq!(limited_thermal.cells(), &limited_cells);
    assert_eq!(limited_thermal.time_s(), 0.0);
    assert_eq!(observer.energy().unwrap().deposited_j, 0.0);
}
