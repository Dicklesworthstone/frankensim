use super::*;

const BACK: CellId = CellId(71);
const BACK_MASS: f64 = 0.017;
const BACK_CONTACT: f64 = 0.13;
const BACK_SOURCE: f64 = 0.3;

fn specimen(
    plateau: bool,
    back_mass: f64,
    back_h: f64,
    back_law: &'static dyn ThermalLaw,
) -> DepositionSimulation<'static> {
    let mut e = events(
        if plateau { [0.0; 2] } else { INITIAL },
        if plateau { [0.0; 2] } else { AUTHORED },
    );
    e.push(TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: BACK,
            segment_id: 5,
            body_id: 91,
            material_id: MaterialId(2),
            mass_kg: back_mass,
            specific_h_j_kg: back_h,
            geometry: CellGeometry {
                start_m: [0.08, 0.0, 0.0],
                end_m: [0.1, 0.0, 0.0],
                volume_m3: AREA * 0.02,
                surface_area_m2: 2.0 * AREA + 4.0 * AREA.sqrt() * 0.02,
            },
            boundary: Boundary {
                ambient_k: 300.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: if plateau { 0.0 } else { BACK_SOURCE },
            },
        }),
    });
    if !plateau {
        for (link, a, b, conductance) in
            [(1, id(0), id(1), CONTACT), (2, id(1), BACK, BACK_CONTACT)]
        {
            e.push(TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: link,
                    a,
                    b,
                    area_m2: AREA,
                    conductance_w_k: conductance,
                }),
            });
        }
    }
    let pore_law: &'static dyn ThermalLaw = if plateau { &PLATEAU } else { &SENSIBLE };
    DepositionSimulation::new(
        e,
        &[
            (MaterialId(0), pore_law),
            (MaterialId(1), pore_law),
            (MaterialId(2), back_law),
        ],
        DepositionBudget {
            max_cells: 16,
            max_links: 16,
            max_events: 32,
            max_substeps: 10000,
        },
    )
    .unwrap()
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0
}

fn harmonic<'a>(
    initial: &DepositionSimulation<'_>,
    materials: &'a [PorousFluidLaw; 2],
) -> TemperaturePorousHeating<'a> {
    let (d, kt, load) = source();
    TemperaturePorousHeating::try_new_oblique_terminated(
        initial,
        &cells(materials),
        gas(),
        d,
        kt,
        load,
        options(),
        || Ok(()),
    )
    .unwrap()
}

// Three independent caloric balances, original rho/K interface fields and
// volume quadrature. Neither native thermal nor load-heating APIs supply this
// reference. Final coordinates integrate backing and the two pore heat rates.
fn reference(materials: &[PorousFluidLaw; 2], steps: usize) -> [f64; 6] {
    let (d, kt, load) = source();
    let rate = |y: [f64; 6]| {
        let p = physical_power([y[0], y[1]], materials, d, kt, load, 64);
        let q01 = CONTACT * (y[1] - y[0]) / CP;
        let qb1 = BACK_CONTACT * (y[2] - y[1]) / CP;
        [
            (p[3] + AUTHORED[0] + q01) / (100.0 * AREA * LENGTHS[0]),
            (p[4] + AUTHORED[1] - q01 + qb1) / (100.0 * AREA * LENGTHS[1]),
            (p[2] + BACK_SOURCE - qb1) / BACK_MASS,
            p[2],
            p[3],
            p[4],
        ]
    };
    let dt = 0.5 / steps as f64;
    let mut y = [INITIAL[0], INITIAL[1], CP * 20.0, 0.0, 0.0, 0.0];
    for _ in 0..steps {
        let a = rate(y);
        let b = rate(std::array::from_fn(|i| y[i] + dt * a[i] * 0.5));
        let c = rate(std::array::from_fn(|i| y[i] + dt * b[i] * 0.5));
        let e = rate(std::array::from_fn(|i| y[i] + dt * c[i]));
        y = std::array::from_fn(|i| y[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + e[i]) / 6.0);
    }
    y
}

#[test]
fn g1_load_heat_and_native_conduction_refine_three_independent_caloric_balances() {
    let materials = materials();
    let initial = specimen(false, BACK_MASS, CP * 20.0, &SENSIBLE);
    let exported = harmonic(&initial, &materials);
    let heating = exported
        .clone()
        .with_load_receiver(&initial, BACK, || Ok(()))
        .unwrap();
    assert_eq!(heating.load_receiver(), Some(BACK));
    assert_eq!(heating.cells(), exported.cells());
    assert_eq!(
        heating.power().load_deposited_power_w,
        heating.power().load_absorbed_power_w
    );
    assert_eq!(heating.power().load_exported_power_w, 0.0);
    assert_eq!(
        exported.power().load_exported_power_w,
        exported.power().load_absorbed_power_w
    );
    let expected = reference(&materials, 512);
    let refined = reference(&materials, 1024);
    for i in 0..6 {
        close(expected[i], refined[i], 3e-9);
    }
    let mut errors = Vec::new();
    for n in [8, 32, 128] {
        let result = heating
            .advance_to(&initial, 0.5, config(0.5 / n as f64), || Ok(()))
            .unwrap();
        let step = &result.step;
        let actual = [
            step.thermal.cells()[&id(0)].specific_h_j_kg,
            step.thermal.cells()[&id(1)].specific_h_j_kg,
            step.thermal.cells()[&BACK].specific_h_j_kg,
            step.acoustic.load_deposited_j,
            step.cell_energy_j[0].1,
            step.cell_energy_j[1].1,
        ];
        errors.push(
            (0..6)
                .map(|i| (actual[i] - expected[i]).abs() / expected[i].abs().max(1.0))
                .fold(0.0, f64::max),
        );
        assert_eq!(
            step.acoustic.load_absorbed_j,
            step.acoustic.load_deposited_j
        );
        assert_eq!(step.acoustic.load_exported_j, 0.0);
        assert_eq!(
            step.load_energy_j,
            Some((BACK, step.acoustic.load_deposited_j))
        );
        let stored: f64 = [id(0), id(1), BACK]
            .iter()
            .map(|id| {
                initial.cells()[id].mass_kg
                    * (step.thermal.cells()[id].specific_h_j_kg
                        - initial.cells()[id].specific_h_j_kg)
            })
            .sum();
        close(
            stored,
            step.acoustic.deposited_j
                + step.acoustic.load_deposited_j
                + 0.5 * (AUTHORED[0] + AUTHORED[1] + BACK_SOURCE),
            2e-11,
        );
        close(step.acoustic.residual_j, 0.0, 3e-10);
        assert_eq!(result.work.gas_evaluations, 2 * result.work.stack_solves);
    }
    assert!(
        errors[1] < 0.3 * errors[0] && errors[2] < 0.3 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 3e-4, "{errors:?}");
    let half = heating
        .advance_to(&initial, 0.25, config(0.03125), || Ok(()))
        .unwrap();
    let split = half
        .heating
        .advance_to(&half.step.thermal, 0.5, config(0.03125), || Ok(()))
        .unwrap();
    let full = heating
        .advance_to(&initial, 0.5, config(0.03125), || Ok(()))
        .unwrap();
    assert_eq!(split.step.thermal.cells(), full.step.thermal.cells());
    close(
        half.step.acoustic.load_deposited_j + split.step.acoustic.load_deposited_j,
        full.step.acoustic.load_deposited_j,
        1e-14,
    );
    let external = exported
        .advance_to(&initial, 0.5, config(0.03125), || Ok(()))
        .unwrap();
    assert_eq!(external.step.acoustic.load_deposited_j, 0.0);
    assert_eq!(
        external.step.acoustic.load_exported_j,
        external.step.acoustic.load_absorbed_j
    );
    assert_eq!(external.step.load_energy_j, None);
    assert!(
        full.step.thermal.temperature_k(BACK).unwrap()
            > external.step.thermal.temperature_k(BACK).unwrap()
    );
    assert_ne!(
        full.heating.power().dissipation,
        external.heating.power().dissipation
    );
}

#[test]
fn g3_frozen_harmonic_and_spectral_backing_heat_share_native_latent_storage() {
    let materials = materials();
    let initial = specimen(true, BACK_MASS, 0.0, &PLATEAU);
    let (d, kt, load) = source();
    let h = harmonic(&initial, &materials)
        .with_load_receiver(&initial, BACK, || Ok(()))
        .unwrap();
    let dt = 0.25;
    let result = h
        .advance_to(&initial, dt, config(dt / 4.0), || Ok(()))
        .unwrap();
    let p = physical_power([0.0; 2], &materials, d, kt, load, 256);
    close(
        result.step.thermal.cells()[&BACK].specific_h_j_kg,
        dt * p[2] / BACK_MASS,
        5e-9,
    );
    assert_eq!(result.step.thermal.temperature_k(BACK).unwrap(), 300.0);
    assert_eq!(result.work.gas_evaluations, 0);
    assert_eq!(result.work.stack_solves, 0);
    let stack = PorousLayerStack::try_new(h.acoustic_layers(), options(), || false).unwrap();
    let excitation = stack
        .prepare_terminated_excitation(d.incident_pressure_pa[0], load, || false)
        .unwrap();
    let partition = [
        [PorousHeatCell {
            cell: id(0),
            end_depth_m: LENGTHS[0],
        }],
        [PorousHeatCell {
            cell: id(1),
            end_depth_m: LENGTHS[1],
        }],
    ];
    let frozen = PorousStackHeating::try_new(
        &excitation,
        &initial,
        &[&partition[0], &partition[1]],
        2,
        || Ok(()),
    )
    .unwrap()
    .with_load_receiver(&initial, BACK, || Ok(()))
    .unwrap();
    assert_eq!(frozen.load_receiver(), Some(BACK));
    let f = frozen
        .advance_to(&initial, dt, config(dt / 4.0), || Ok(()))
        .unwrap();
    assert_eq!(f.thermal.cells(), result.step.thermal.cells());
    assert_eq!(f.acoustic, result.step.acoustic);
    let mut d2 = d;
    d2.angular_frequency_rad_s *= 1.4;
    d2.incident_pressure_pa[0] = d2.incident_pressure_pa[0].scale(0.6);
    let load2 = PassivePorousTermination::try_new_impedance(
        C64::new(0.4, 0.25).scale(d.reference_impedance_pa_s_m3),
    )
    .unwrap();
    for bins in 1..=2 {
        let drives = [d, d2];
        let kt_values = [kt, 0.0];
        let loads = [load, load2];
        let spectral = PorousSpectrumHeating::try_new_oblique_terminated(
            &initial,
            &cells(&materials),
            gas(),
            &drives[..bins],
            &kt_values[..bins],
            &loads[..bins],
            PorousSpectrumOptions::default(),
            || Ok(()),
        )
        .unwrap()
        .with_load_receiver(&initial, BACK, || Ok(()))
        .unwrap();
        let covariance: Vec<_> = drives[..bins]
            .iter()
            .map(|d| PorousCovarianceBin {
                angular_frequency_rad_s: d.angular_frequency_rad_s,
                area_m2: AREA,
                reference_impedance_pa_s_m3: d.reference_impedance_pa_s_m3,
                covariance: PressureCovariance::from_peak_pressure(d.incident_pressure_pa).unwrap(),
            })
            .collect();
        let modes: Vec<_> = covariance
            .iter()
            .enumerate()
            .map(|(i, &bin)| PorousAngularMode {
                bin,
                tangential_wavenumber_rad_m: kt_values[i],
            })
            .collect();
        let cov = PorousCrossSpectrumHeating::try_new_oblique_terminated(
            &initial,
            &cells(&materials),
            gas(),
            &covariance,
            &kt_values[..bins],
            &loads[..bins],
            PorousSpectrumOptions::default(),
            || Ok(()),
        )
        .unwrap()
        .with_load_receiver(&initial, BACK, || Ok(()))
        .unwrap();
        let angular = PorousCrossSpectrumHeating::try_new_independent_modes_terminated(
            &initial,
            &cells(&materials),
            gas(),
            &modes,
            &loads[..bins],
            PorousSpectrumOptions::default(),
            || Ok(()),
        )
        .unwrap()
        .with_load_receiver(&initial, BACK, || Ok(()))
        .unwrap();
        assert_eq!(spectral.load_receiver(), Some(BACK));
        assert_eq!(cov.load_receiver(), Some(BACK));
        assert_eq!(angular.load_receiver(), Some(BACK));
        let s = spectral
            .advance_to(&initial, dt, config(dt / 4.0), || Ok(()))
            .unwrap();
        let c = cov
            .advance_to(&initial, dt, config(dt / 4.0), || Ok(()))
            .unwrap();
        let a = angular
            .advance_to(&initial, dt, config(dt / 4.0), || Ok(()))
            .unwrap();
        assert_eq!(s.step.thermal.cells(), c.step.thermal.cells());
        assert_eq!(s.step.thermal.cells(), a.step.thermal.cells());
        assert_eq!(s.step.acoustic, c.step.acoustic);
        assert_eq!(s.step.acoustic, a.step.acoustic);
        if bins == 1 {
            assert_eq!(s.step.thermal.cells(), result.step.thermal.cells());
        }
        let expected: f64 = (0..bins)
            .map(|i| {
                physical_power([0.0; 2], &materials, drives[i], kt_values[i], loads[i], 256)[2] * dt
            })
            .sum();
        close(s.step.acoustic.load_deposited_j, expected, 4e-10);
        close(
            s.bin_energy.iter().map(|e| e.load_deposited_j).sum(),
            s.step.acoustic.load_deposited_j,
            1e-14,
        );
        assert_eq!(s.step.acoustic.load_exported_j, 0.0);
        close(
            s.step.thermal.cells()[&BACK].specific_h_j_kg * BACK_MASS,
            expected,
            4e-10,
        );
    }
}

#[test]
fn g0_g4_backing_receiver_refuses_alias_rebinding_and_domain_failure_without_publication() {
    let materials = materials();
    let initial = specimen(true, BACK_MASS, 0.0, &PLATEAU);
    let h = harmonic(&initial, &materials);
    assert!(
        h.clone()
            .with_load_receiver(&initial, id(0), || Ok(()))
            .is_err()
    );
    assert!(
        h.clone()
            .with_load_receiver(&initial, CellId(999), || Ok(()))
            .is_err()
    );
    assert!(matches!(
        h.clone()
            .with_load_receiver(&initial, BACK, || Err(DepositionError::Cancelled)),
        Err(PorousHeatingError::Thermal(DepositionError::Cancelled))
    ));
    let bound = h.with_load_receiver(&initial, BACK, || Ok(())).unwrap();
    assert!(
        bound
            .clone()
            .with_load_receiver(&initial, BACK, || Ok(()))
            .is_err()
    );
    let foreign = specimen(true, BACK_MASS * 2.0, 0.0, &PLATEAU);
    assert!(
        bound
            .advance_to(&foreign, 0.1, config(0.01), || Ok(()))
            .is_err()
    );
    let exhausted = specimen(true, BACK_MASS, 1e4 - 1e-4, &PLATEAU);
    assert!(matches!(
        bound.advance_to(&exhausted, 0.1, config(0.1), || Ok(())),
        Err(PorousHeatingError::Thermal(DepositionError::MaterialDomain))
    ));
    assert_eq!(exhausted.cells()[&BACK].specific_h_j_kg, 1e4 - 1e-4);
    let before = initial.cells().clone();
    let mut checks = 0;
    assert!(
        bound
            .advance_to(&initial, 0.1, config(0.01), || {
                checks += 1;
                if checks == 30 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    assert_eq!(initial.cells(), &before);
    let (d, kt, _) = source();
    let rigid = TemperaturePorousHeating::try_new_oblique_terminated(
        &initial,
        &cells(&materials),
        gas(),
        d,
        kt,
        PassivePorousTermination::rigid(),
        options(),
        || Ok(()),
    )
    .unwrap()
    .with_load_receiver(&initial, BACK, || Ok(()))
    .unwrap();
    let step = rigid
        .advance_to(&initial, 0.1, config(0.025), || Ok(()))
        .unwrap();
    assert_eq!(step.step.acoustic.load_absorbed_j, 0.0);
    assert_eq!(step.step.load_energy_j, Some((BACK, 0.0)));
    assert_eq!(step.step.thermal.cells()[&BACK].specific_h_j_kg, 0.0);
}
