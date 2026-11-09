//! Original volume injection, analytic compression and transactional refresh.

use super::*;
use fs_feec::acoustics::helmholtz::AcousticFeError;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

fn injected(drive: InteriorTemperatureVolumeDrive) -> InteriorTemperatureAcoustics {
    let mut value = problem(false);
    value.inward_velocity.clear();
    value.volume_drives.push(drive);
    value
}

#[test]
fn uniform_original_injection_has_exact_compression_work_and_plateau_cache_reuse() {
    let (mesh, positions) = fixture();
    let geometry = HelmholtzFe::new(
        &mesh,
        &positions,
        2,
        Medium {
            density_kg_m3: 1.0,
            sound_speed_m_s: 2.0,
        },
    )
    .unwrap();
    let law = Plateau;
    let porous = laws();
    let mut cells = mapping(&law, Some(&porous));
    cells[1].porous = Some(&porous[0]);
    let thermal = initial(events([0.0; 2], [0.0; 2], false), &law);
    let volume = (VOLUMES[0] + VOLUMES[1]) * SCALE.powi(3);
    let calls = Arc::new(AtomicUsize::new(0));
    let sampled = Arc::clone(&calls);
    let q = C64::new(1e-5, -3e-6);
    let drive = InteriorTemperatureVolumeDrive {
        density_per_m3: Arc::new(move |_, _| {
            sampled.fetch_add(1, Ordering::Relaxed);
            Ok(1.0 / volume)
        }),
        quadrature_order: 4,
        amplitude_m3_s: q,
    };
    let cloned = drive.clone();
    assert!(Arc::ptr_eq(&drive.density_per_m3, &cloned.density_per_m3));
    assert!(format!("{drive:?}").contains("InteriorTemperatureVolumeDrive"));
    assert_eq!(calls.load(Ordering::Relaxed), 0);
    let heating = TemperatureInteriorHeating::try_new(
        &geometry,
        &thermal,
        &cells,
        &[],
        fluid([280.0, 500.0]),
        injected(drive),
        2,
        || Ok(()),
    )
    .unwrap();
    let samples = calls.load(Ordering::Relaxed);
    assert!(samples > 0);

    // Rigid homogeneous cavity: i*omega*(V/K)*p = Q on exp(+i*omega*t).
    // This independent compression oracle retains the original porous K;
    // pressure is constant in all vertex coordinates and zero on P2 edges.
    let gas = GasState::try_new(&SPEC, 300.0, 101325.0).unwrap();
    let omega = core::f64::consts::TAU * FREQUENCY;
    let bulk = porous[0]
        .effective_fluid(&gas, omega)
        .unwrap()
        .dynamic_bulk_modulus_pa();
    let product = bulk * q;
    let pressure = C64::new(product.im, -product.re).scale(1.0 / (omega * volume));
    let watts = 0.5 * bulk.im * q.norm_sq() / (omega * volume);
    assert!(watts > 0.0);
    for (i, &actual) in heating.field().pressure_pa.iter().enumerate() {
        let expected = if i < positions.len() {
            pressure
        } else {
            C64::ZERO
        };
        assert!((actual - expected).abs() < 2e-8 * pressure.abs());
    }
    close(heating.field().volume_source_power_w, watts, 2e-8 * watts);
    close(
        heating.field().compressional_dissipation_w,
        watts,
        2e-8 * watts,
    );
    let result = heating
        .advance_to(&thermal, 2.0 * DT, config(DT), || Ok(()))
        .unwrap();
    assert_eq!(result.work, InteriorFeedbackWork::default());
    assert_eq!(calls.load(Ordering::Relaxed), samples);
    assert!(core::ptr::eq(heating.prepared(), result.heating.prepared()));
    assert!(core::ptr::eq(heating.field(), result.heating.field()));
    close(
        result.step.acoustic.compressional_dissipation_j,
        watts * 2.0 * DT,
        3e-8 * watts * 2.0 * DT,
    );
    close(
        result.step.acoustic.deposited_j,
        watts * 2.0 * DT,
        3e-8 * watts * 2.0 * DT,
    );
    close(
        result.step.heat_report.external_heat_j,
        result.step.acoustic.deposited_j,
        2e-12 * watts * 2.0 * DT,
    );
    for (i, &(id, energy)) in result.step.cell_energy_j.iter().enumerate() {
        assert_eq!(id, IDS[i]);
        close(
            energy,
            watts * 2.0 * DT * VOLUMES[i] / VOLUMES.iter().sum::<f64>(),
            3e-8 * watts * 2.0 * DT,
        );
        assert_eq!(result.step.thermal.temperature_k(id).unwrap(), 300.0);
        assert!(result.step.thermal.cells()[&id].specific_h_j_kg > 0.0);
    }
}

#[test]
fn volume_admission_and_callback_refusal_preserve_original_feedback_for_retry() {
    let (mesh, positions) = fixture();
    let geometry = geometry(&mesh, &positions);
    let law = caloric();
    let porous = laws();
    let thermal = initial(events([0.0, 40.0 * CP], [0.12, 0.1], false), &law);
    let calls = Arc::new(AtomicUsize::new(0));
    let unavailable = Arc::new(AtomicBool::new(false));
    let sampled = Arc::clone(&calls);
    let fail = Arc::clone(&unavailable);
    let drive = InteriorTemperatureVolumeDrive {
        density_per_m3: Arc::new(move |element, position| {
            sampled.fetch_add(1, Ordering::Relaxed);
            if fail.load(Ordering::Relaxed) {
                return Err(HeterogeneousError::Acoustic(
                    AcousticFeError::InvalidInput {
                        what: "original volume profile unavailable",
                    },
                ));
            }
            Ok((1.0 + 0.2 * element as f64 + position[0] / SCALE) / (0.5 * SCALE.powi(3)))
        }),
        quadrature_order: 4,
        amplitude_m3_s: C64::new(1e-5, -4e-6),
    };
    let source = fluid([280.0, 500.0]);
    let material_calls = source.calls.clone();
    let mut excessive = injected(drive.clone());
    excessive.volume_drives = vec![drive.clone(); MAX_VOLUME_SOURCES + 1];
    assert!(matches!(TemperatureInteriorHeating::try_new(
        &geometry, &thermal, &mapping(&law, Some(&porous)), &[], source.clone(),
        excessive, 2, || Ok(()),
    ), Err(InteriorFeedbackError::Heating(InteriorHeatingError::Acoustic(
        HeterogeneousError::Limit { requested, maximum, .. }
    ))) if requested == MAX_VOLUME_SOURCES + 1 && maximum == MAX_VOLUME_SOURCES));
    for amplitude in [C64::new(f64::NAN, 0.0), C64::new(0.0, f64::INFINITY)] {
        let mut bad = drive.clone();
        bad.amplitude_m3_s = amplitude;
        assert!(matches!(
            TemperatureInteriorHeating::try_new(
                &geometry,
                &thermal,
                &mapping(&law, Some(&porous)),
                &[],
                source.clone(),
                injected(bad),
                2,
                || Ok(()),
            ),
            Err(InteriorFeedbackError::Heating(
                InteriorHeatingError::Invalid(_)
            ))
        ));
    }
    assert_eq!(material_calls.get(), 0);
    assert_eq!(calls.load(Ordering::Relaxed), 0);
    for bad in [
        InteriorTemperatureVolumeDrive {
            quadrature_order: 2,
            ..drive.clone()
        },
        InteriorTemperatureVolumeDrive {
            density_per_m3: Arc::new(|_, _| Ok(f64::NAN)),
            ..drive.clone()
        },
    ] {
        assert!(matches!(
            TemperatureInteriorHeating::try_new(
                &geometry,
                &thermal,
                &mapping(&law, Some(&porous)),
                &[],
                source.clone(),
                injected(bad),
                2,
                || Ok(()),
            ),
            Err(InteriorFeedbackError::Heating(
                InteriorHeatingError::Acoustic(_)
            ))
        ));
    }
    assert_eq!(calls.load(Ordering::Relaxed), 0);
    let heating = TemperatureInteriorHeating::try_new(
        &geometry,
        &thermal,
        &mapping(&law, Some(&porous)),
        &[],
        source,
        injected(drive),
        2,
        || Ok(()),
    )
    .unwrap();
    let samples = calls.load(Ordering::Relaxed);
    let original = thermal.cells().clone();
    let pressure = heating.field().pressure_pa.clone();
    let accepted = heating
        .advance_to(&thermal, 2.0 * DT, config(DT), || Ok(()))
        .unwrap();
    assert_eq!(accepted.work.field_solves, 2);
    assert_eq!(accepted.work.fluid_evaluations, 4);
    assert_eq!(calls.load(Ordering::Relaxed), 3 * samples);
    assert_ne!(accepted.heating.field().pressure_pa, pressure);

    // A profile can refuse an endpoint rebuild; its values have not changed.
    // The accepted thermal cells and old field must remain a reusable pair.
    unavailable.store(true, Ordering::Relaxed);
    assert!(matches!(
        heating.advance_to(&thermal, 2.0 * DT, config(DT), || Ok(())),
        Err(InteriorFeedbackError::Heating(
            InteriorHeatingError::Acoustic(HeterogeneousError::Acoustic(
                AcousticFeError::InvalidInput {
                    what: "original volume profile unavailable",
                }
            ))
        ))
    ));
    unavailable.store(false, Ordering::Relaxed);
    assert_eq!(thermal.cells(), &original);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(heating.field().pressure_pa, pressure);

    // Cancel only after the next source projection has actually sampled.
    let before = calls.load(Ordering::Relaxed);
    assert!(matches!(
        heating.advance_to(&thermal, 2.0 * DT, config(DT), || {
            if calls.load(Ordering::Relaxed) > before {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(InteriorFeedbackError::Heating(
            InteriorHeatingError::Thermal(DepositionError::Cancelled)
        ))
    ));
    assert_eq!(thermal.cells(), &original);
    assert_eq!(heating.field().pressure_pa, pressure);
    let retry = heating
        .advance_to(&thermal, 2.0 * DT, config(DT), || Ok(()))
        .unwrap();
    assert_eq!(retry.step.thermal.cells(), accepted.step.thermal.cells());
    assert_eq!(
        retry.heating.field().pressure_pa,
        accepted.heating.field().pressure_pa
    );
    assert_eq!(retry.step.acoustic, accepted.step.acoustic);
    assert_eq!(retry.work, accepted.work);
}
