//! Coherent original boundary/wall drives across accepted-temperature refresh.
//! References prepare the original material/FE owners and apply their positive
//! loss forms through the existing frozen native thermal interface.

use super::*;
use fs_feec::acoustics::helmholtz::heterogeneous::time::mechanical::traction::WallTractionProfile;
use fs_feec::acoustics::helmholtz::heterogeneous::{
    EquivalentFluidPressureField, MAX_BOUNDARY_SOURCES, MAX_VOLUME_SOURCES,
};

mod volume;

const DT: f64 = 0.03125;

fn driven_problem() -> InteriorTemperatureAcoustics {
    let mut value = problem(true);
    value.inward_velocity.clear();
    value.boundary_drives.push(InteriorTemperatureBoundaryDrive {
        // A true midside-only trace; its three corner values are all zero.
        traces: vec![([0, 2, 3], [0.0, 0.0, 0.0, 1.0, 0.0, 0.0])],
        quadrature_order: 4,
        amplitude_m_s: C64::new(0.4, -0.1),
    });
    value.traction_drives.push(InteriorTemperatureTractionDrive {
        profile: WallTractionProfile {
            wall: 0,
            coefficients: vec![(0, 1.0), (1, 1.0), (3, 1.0)],
        },
        amplitude_pa: C64::new(700.0, -150.0),
    });
    value
}

fn native_field(
    prepared: &PreparedEquivalentFluid<'_, '_>,
    problem: &InteriorTemperatureAcoustics,
) -> EquivalentFluidPressureField {
    let boundary: Vec<_> = problem
        .boundary_drives
        .iter()
        .map(|drive| {
            prepared
                .quadratic_velocity_source_with_quadrature(
                    &drive.traces,
                    drive.quadrature_order,
                    &mut || false,
                )
                .unwrap()
        })
        .collect();
    let traction: Vec<_> = problem
        .traction_drives
        .iter()
        .map(|drive| prepared.prepare_wall_traction(&drive.profile, &mut || false).unwrap())
        .collect();
    let volume: Vec<_> = problem
        .volume_drives
        .iter()
        .map(|drive| {
            prepared
                .volume_velocity_source(
                    drive.quadrature_order,
                    &mut |element, position| (drive.density_per_m3)(element, position),
                    &mut || false,
                )
                .unwrap()
        })
        .collect();
    let boundary_values: Vec<_> = boundary
        .iter()
        .zip(&problem.boundary_drives)
        .map(|(source, drive)| (source, drive.amplitude_m_s))
        .collect();
    let traction_values: Vec<_> = traction
        .iter()
        .zip(&problem.traction_drives)
        .map(|(source, drive)| (source, drive.amplitude_pa))
        .collect();
    let volume_values: Vec<_> = volume
        .iter()
        .zip(&problem.volume_drives)
        .map(|(source, drive)| (source, drive.amplitude_m3_s))
        .collect();
    prepared
        .solve_with_sources(
            &problem.inward_velocity,
            &boundary_values,
            &volume_values,
            &traction_values,
            &mut || false,
        )
        .unwrap()
}

fn powers(field: &EquivalentFluidPressureField) -> [f64; 4] {
    [
        field.density_dissipation_w,
        field.compressional_dissipation_w,
        field.wall_accounting.resistance_dissipation_w,
        field.wall_accounting.relaxation_dissipation_w,
    ]
}

#[test]
fn mixed_quadratic_volume_and_traction_drives_refresh_on_original_materials_and_heat_once() {
    let (mesh, positions) = fixture();
    let geometry = HelmholtzFe::new(
        &mesh,
        &positions,
        2,
        Medium { density_kg_m3: 1.0, sound_speed_m_s: 2.0 },
    )
    .unwrap();
    let law = caloric();
    let porous = laws();
    let thermal = initial(events([0.0, 40.0 * CP], [0.12, 0.1], true), &law);
    let original_cells = thermal.cells().clone();
    let source = fluid([280.0, 500.0]);
    let calls = source.calls.clone();
    let mut problem = driven_problem();
    problem.volume_drives.push(InteriorTemperatureVolumeDrive {
        density_per_m3: std::sync::Arc::new(|element, position| {
            let sign = if element == 0 { 1.0 } else { -0.25 };
            Ok(sign * (1.0 + position[0] / SCALE) / (0.5 * SCALE.powi(3)))
        }),
        quadrature_order: 4,
        amplitude_m3_s: C64::new(1e-5, -4e-6),
    });
    let heating = TemperatureInteriorHeating::try_new(
        &geometry, &thermal, &mapping(&law, Some(&porous)), &wall_routes(&law),
        source, problem.clone(), 3, || Ok(()),
    )
    .unwrap();
    assert_eq!(calls.get(), 2);
    assert!(heating.field().pressure_pa.iter().any(|p| p.norm_sq() > 0.0));
    assert!(powers(heating.field()).iter().all(|&power| power > 0.0));

    let result = heating.advance_to(&thermal, 2.0 * DT, config(DT), || Ok(())).unwrap();
    assert_eq!(result.step.heat_report.substeps, 2);
    assert_eq!(result.work.field_solves, 2);
    assert_eq!(result.work.fluid_evaluations, 4);
    assert_eq!(result.work.porous_evaluations, 4);
    assert_eq!(calls.get(), 6);

    // Rebuild the ORIGINAL source laws at each native accepted temperature.
    // Hold each independently prepared field only over its following interval;
    // neither endpoint watts nor signed source work are an extra heat source.
    let mut reference = thermal.clone();
    let mut expected = [0.0; 4];
    let mut cell_j = [0.0; 2];
    for index in 0..=2 {
        let states = IDS.map(|id| {
            GasState::try_new(&SPEC, reference.temperature_k(id).unwrap(), 101325.0).unwrap()
        });
        let regions = core::array::from_fn::<_, 2, _>(|i| AcousticVolumeRegion::Porous {
            law: &porous[i],
            fluid: &states[i],
        });
        let original = PorousRegionsInterior::new_with_coupled_mechanical_walls(
            &geometry, FREQUENCY, &regions, &[0, 1], &problem.admittances,
            &problem.mechanical_walls, &problem.pressure_release_faces,
            problem.controls, &mut || false,
        )
        .unwrap();
        let field = native_field(original.operator(), &problem);
        assert!(field.power_balance_residual_w().abs() < 1e-9);
        if index == 0 {
            assert_eq!(field.pressure_pa, heating.field().pressure_pa);
            assert_eq!(
                field.walls[0].normal_velocity_m_s,
                heating.field().walls[0].normal_velocity_m_s,
            );
        }
        if index == 2 {
            for (actual, expected) in result.heating.field().pressure_pa.iter().zip(&field.pressure_pa) {
                assert!((*actual - *expected).abs() < 1e-8);
            }
            for (actual, expected) in result.heating.field().walls[0].normal_velocity_m_s.iter()
                .zip(&field.walls[0].normal_velocity_m_s)
            {
                assert!((*actual - *expected).abs() < 1e-10);
            }
            break;
        }
        for (energy, watts) in expected.iter_mut().zip(powers(&field)) {
            *energy += DT * watts;
        }
        let held = InteriorAcousticHeating::try_new_with_mechanical_walls(
            original.operator(), &field.pressure_pa, &field.walls, &reference,
            &routes(&law), &wall_routes(&law), 3, || Ok(()),
        )
        .unwrap();
        let step = held.advance_to(&reference, (index + 1) as f64 * DT, config(DT), || Ok(())).unwrap();
        for (i, &(id, joules)) in step.cell_energy_j.iter().enumerate() {
            assert_eq!(id, IDS[i]);
            cell_j[i] += joules;
        }
        reference = step.thermal;
    }
    let actual = result.step.acoustic;
    for (actual, expected) in [
        actual.density_dissipation_j,
        actual.compressional_dissipation_j,
        actual.wall_resistance_dissipation_j,
        actual.wall_relaxation_dissipation_j,
    ].into_iter().zip(expected) {
        close(actual, expected, 2e-12);
    }
    close(actual.deposited_j, expected.iter().sum(), 2e-12);
    for (i, &(id, joules)) in result.step.cell_energy_j.iter().enumerate() {
        assert_eq!(id, IDS[i]);
        close(joules, cell_j[i], 2e-12);
    }
    for id in IDS.into_iter().chain([WALL]) {
        close(result.step.thermal.cells()[&id].specific_h_j_kg,
            reference.cells()[&id].specific_h_j_kg, 2e-10);
        assert_eq!(result.step.thermal.cells()[&id].mass_kg, original_cells[&id].mass_kg);
        assert_eq!(result.step.thermal.cells()[&id].geometry, original_cells[&id].geometry);
    }
    close(result.step.heat_report.external_heat_j,
        actual.deposited_j + (0.12 + 0.1 + 0.05) * 2.0 * DT, 2e-11);
    close(result.step.heat_report.residual_j, 0.0, 2e-11);
    assert_ne!(result.heating.field().pressure_pa, heating.field().pressure_pa);
    for (i, state) in result.heating.fluid_states().iter().enumerate() {
        assert_eq!(state.temperature, result.step.thermal.temperature_k(IDS[i]).unwrap());
    }
    assert_eq!(thermal.cells(), &original_cells);
}

#[test]
fn mechanical_source_alone_heats_a_plateau_without_rebuilding_its_native_owner() {
    let (mesh, positions) = fixture();
    let geometry = geometry(&mesh, &positions);
    let law = Plateau;
    let porous = laws();
    let thermal = initial(events([0.0; 2], [0.0; 2], true), &law);
    let source = fluid([280.0, 500.0]);
    let calls = source.calls.clone();
    let mut problem = driven_problem();
    problem.boundary_drives.clear();
    assert!(problem.inward_velocity.is_empty());
    let heating = TemperatureInteriorHeating::try_new(
        &geometry, &thermal, &mapping(&law, Some(&porous)), &wall_routes(&law),
        source, problem, 3, || Ok(()),
    )
    .unwrap();
    assert!(heating.field().pressure_pa.iter().any(|p| p.norm_sq() > 0.0));
    assert!(powers(heating.field()).iter().all(|&power| power > 0.0));
    let held = InteriorAcousticHeating::try_new_with_mechanical_walls(
        heating.prepared(), &heating.field().pressure_pa, &heating.field().walls,
        &thermal, &routes(&law), &wall_routes(&law), 3, || Ok(()),
    )
    .unwrap();
    let expected = held.advance_to(&thermal, 2.0 * DT, config(DT), || Ok(())).unwrap();
    let result = heating.advance_to(&thermal, 2.0 * DT, config(DT), || Ok(())).unwrap();
    assert_eq!(result.work, InteriorFeedbackWork::default());
    assert_eq!(calls.get(), 2);
    assert!(core::ptr::eq(heating.prepared(), result.heating.prepared()));
    assert!(core::ptr::eq(heating.field(), result.heating.field()));
    assert_eq!(result.step.thermal.cells(), expected.thermal.cells());
    assert_eq!(result.step.acoustic, expected.acoustic);
    assert_eq!(result.step.cell_energy_j, expected.cell_energy_j);
    close(result.step.acoustic.deposited_j,
        powers(heating.field()).iter().sum::<f64>() * 2.0 * DT, 2e-12);
    close(result.step.heat_report.external_heat_j,
        result.step.acoustic.deposited_j + 0.05 * 2.0 * DT, 2e-11);
    for id in IDS.into_iter().chain([WALL]) {
        assert_eq!(result.step.thermal.temperature_k(id).unwrap(), 300.0);
        assert!(result.step.thermal.cells()[&id].specific_h_j_kg > 0.0);
    }
}

#[test]
fn source_admission_and_late_feedback_cancellation_preserve_the_original_snapshot() {
    let (mesh, positions) = fixture();
    let geometry = geometry(&mesh, &positions);
    let law = caloric();
    let porous = laws();
    let thermal = initial(events([0.0, 40.0 * CP], [0.12, 0.1], true), &law);
    let source = fluid([280.0, 500.0]);
    let calls = source.calls.clone();
    for traction in [false, true] {
        let mut problem = driven_problem();
        let maximum = if traction { MAX_VOLUME_SOURCES } else { MAX_BOUNDARY_SOURCES };
        if traction {
            problem.traction_drives = vec![problem.traction_drives[0].clone(); maximum + 1];
        } else {
            problem.boundary_drives = vec![problem.boundary_drives[0].clone(); maximum + 1];
        }
        assert!(matches!(
            TemperatureInteriorHeating::try_new(
                &geometry, &thermal, &mapping(&law, Some(&porous)), &wall_routes(&law),
                source.clone(), problem, 3, || Ok(()),
            ),
            Err(InteriorFeedbackError::Heating(InteriorHeatingError::Acoustic(
                HeterogeneousError::Limit { requested, maximum: cap, .. }
            ))) if requested == maximum + 1 && cap == maximum
        ));
        assert_eq!(calls.get(), 0);
    }
    for traction in [false, true] {
        let mut problem = driven_problem();
        if traction {
            problem.traction_drives[0].amplitude_pa = C64::new(f64::NAN, 0.0);
        } else {
            problem.boundary_drives[0].amplitude_m_s = C64::new(0.0, f64::INFINITY);
        }
        assert!(matches!(
            TemperatureInteriorHeating::try_new(
                &geometry, &thermal, &mapping(&law, Some(&porous)), &wall_routes(&law),
                source.clone(), problem, 3, || Ok(()),
            ),
            Err(InteriorFeedbackError::Heating(InteriorHeatingError::Invalid(_)))
        ));
        assert_eq!(calls.get(), 0);
    }
    for bad_wall_index in [false, true] {
        let mut problem = driven_problem();
        if bad_wall_index {
            problem.traction_drives[0].profile.wall = 1;
        } else {
            // Native quadratic inward motion cannot share a pressure-release face.
            problem.pressure_release_faces.push([0, 2, 3]);
        }
        assert!(matches!(
            TemperatureInteriorHeating::try_new(
                &geometry, &thermal, &mapping(&law, Some(&porous)), &wall_routes(&law),
                source.clone(), problem, 3, || Ok(()),
            ),
            Err(InteriorFeedbackError::Heating(InteriorHeatingError::Acoustic(_)))
        ));
    }
    let heating = TemperatureInteriorHeating::try_new(
        &geometry, &thermal, &mapping(&law, Some(&porous)), &wall_routes(&law),
        source, driven_problem(), 3, || Ok(()),
    )
    .unwrap();
    let original_cells = thermal.cells().clone();
    let pressure = heating.field().pressure_pa.clone();
    let mut polls = 0;
    let accepted = heating.advance_to(&thermal, 2.0 * DT, config(DT), || {
        polls += 1;
        Ok(())
    }).unwrap();
    let mut cancelled_polls = 0;
    assert!(matches!(
        heating.advance_to(&thermal, 2.0 * DT, config(DT), || {
            cancelled_polls += 1;
            if cancelled_polls >= polls * 3 / 4 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(InteriorFeedbackError::Heating(
            InteriorHeatingError::Thermal(DepositionError::Cancelled)
        ))
    ));
    assert_eq!(thermal.cells(), &original_cells);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(heating.field().pressure_pa, pressure);
    let retry = heating.advance_to(&thermal, 2.0 * DT, config(DT), || Ok(())).unwrap();
    assert_eq!(retry.step.thermal.cells(), accepted.step.thermal.cells());
    assert_eq!(retry.heating.field().pressure_pa, accepted.heating.field().pressure_pa);
    assert_eq!(retry.step.acoustic, accepted.step.acoustic);
    assert_eq!(retry.work, accepted.work);
}
