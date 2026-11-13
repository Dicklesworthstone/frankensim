//! Native covariance heat on explicit numerical caloric fixtures. These checks
//! exercise model composition, not measured qualification of a real specimen.
use fs_couple::porous_heating::PorousHeatingError;
use fs_couple::porous_heating::feedback::covariance::*;
use fs_couple::porous_heating::feedback::spectrum::*;
use fs_couple::porous_heating::feedback::*;
use fs_duct::porous::PressureCovariance;
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_math::c64::C64;
use fs_thermal::enthalpy::*;

#[path = "porous_feedback_covariance/angular.rs"]
mod angular;
#[path = "porous_feedback_covariance/refinement.rs"]
mod refinement;

const AREA: f64 = 0.00390625;
const LENGTHS: [f64; 2] = [0.03125, 0.046875];
const IDS: [CellId; 2] = [CellId(37), CellId(9)];
const CP: f64 = 200.0;
static SPEC: GasSpec = GasSpec::dry_air_ussa1976();
struct Caloric(bool);
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if self.0 {
            if (0.0..=10000.0).contains(&h) {
                Ok(300.0)
            } else {
                Err(DepositionError::MaterialDomain)
            }
        } else {
            let t = 300.0 + h / CP;
            if (280.0..=500.0).contains(&t) {
                Ok(t)
            } else {
                Err(DepositionError::MaterialDomain)
            }
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        if self.0 { 0.0 } else { 1.0 / CP }
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([80 + u8::from(self.0); 32])
    }
}
static SENSIBLE: Caloric = Caloric(false);
static PLATEAU: Caloric = Caloric(true);
fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-9,
    }
}
fn options() -> PorousSpectrumOptions {
    PorousSpectrumOptions {
        maximum_bins: 4,
        maximum_cell_bins: 16,
        ..PorousSpectrumOptions::default()
    }
}
fn gas() -> IsobaricPorousGas<'static> {
    IsobaricPorousGas::try_new(&SPEC, 101325.0, [280.0, 450.0]).unwrap()
}
fn reference() -> f64 {
    GasState::try_new(&SPEC, 300.0, 101325.0)
        .unwrap()
        .characteristic_impedance
        / AREA
}
fn materials() -> [PorousFluidLaw; 2] {
    std::array::from_fn(|i| {
        let mut card = PorousCard::melamine_uf_nguyen2024();
        card.name = "authored covariance absorption fixture";
        card.provenance = "numerical correlated-noise test, not measured qualification";
        card.flow_resistivity *= [1.0, 1.7][i];
        PorousFluidLaw::try_new(card, 1.81e-5).unwrap()
    })
}
fn routes(materials: &[PorousFluidLaw; 2]) -> [PorousTemperatureCell<'_>; 2] {
    std::array::from_fn(|i| PorousTemperatureCell {
        cell: IDS[i],
        material: &materials[i],
    })
}
fn covariance(a: f64, b: f64, c: C64) -> PressureCovariance {
    PressureCovariance::try_new([C64::from_re(a), c, c.conj(), C64::from_re(b)]).unwrap()
}
fn bins() -> [PorousCovarianceBin; 2] {
    [
        (
            1100.0,
            covariance(15000.0, 22000.0, C64::new(-4000.0, 6000.0)),
        ),
        (
            450.0,
            covariance(30000.0, 16000.0, C64::new(5000.0, 8000.0)),
        ),
    ]
    .map(|(frequency, covariance)| PorousCovarianceBin {
        angular_frequency_rad_s: core::f64::consts::TAU * frequency,
        area_m2: AREA,
        reference_impedance_pa_s_m3: reference(),
        covariance,
    })
}
fn events() -> Vec<TimedEvent> {
    let mut start = 0.0;
    let mut result = Vec::new();
    for i in 0..2 {
        let end = start + LENGTHS[i];
        result.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
                segment_id: i as u64,
                body_id: 81,
                material_id: MaterialId(i as u64),
                mass_kg: 100.0 * AREA * LENGTHS[i],
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [start, 0.0, 0.0],
                    end_m: [end, 0.0, 0.0],
                    volume_m3: AREA * LENGTHS[i],
                    surface_area_m2: 2.0 * AREA + 4.0 * AREA.sqrt() * LENGTHS[i],
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: [0.2, -0.05][i],
                },
            }),
        });
        start = end;
    }
    result.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 1,
            a: IDS[0],
            b: IDS[1],
            area_m2: AREA,
            conductance_w_k: 0.02,
        }),
    });
    result
}
fn thermal(
    events: Vec<TimedEvent>,
    laws: [&'static dyn ThermalLaw; 2],
) -> DepositionSimulation<'static> {
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), laws[0]), (MaterialId(1), laws[1])],
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 16,
            max_substeps: 8192,
        },
    )
    .unwrap()
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0
}
fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= 2e-14 + relative * a.abs().max(b.abs()),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn g1_covariance_preserves_coherent_spectrum_exactly_and_adds_uncorrelated_face_powers() {
    let materials = materials();
    let initial = thermal(events(), [&SENSIBLE; 2]);
    let p = [
        [C64::new(150.0, -30.0), C64::new(40.0, 60.0)],
        [C64::new(80.0, 100.0), C64::new(-20.0, 30.0)],
    ];
    let phasors: Vec<_> = bins()
        .iter()
        .zip(p)
        .map(|(b, p)| PorousHarmonicDrive {
            angular_frequency_rad_s: b.angular_frequency_rad_s,
            area_m2: b.area_m2,
            reference_impedance_pa_s_m3: b.reference_impedance_pa_s_m3,
            incident_pressure_pa: p,
        })
        .collect();
    let coherent: Vec<_> = bins()
        .iter()
        .zip(p)
        .map(|(b, p)| PorousCovarianceBin {
            covariance: PressureCovariance::from_peak_pressure(p).unwrap(),
            ..*b
        })
        .collect();
    let old = PorousSpectrumHeating::try_new(
        &initial,
        &routes(&materials),
        gas(),
        &phasors,
        options(),
        || Ok(()),
    )
    .unwrap();
    let new = PorousCrossSpectrumHeating::try_new(
        &initial,
        &routes(&materials),
        gas(),
        &coherent,
        options(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(old.power(), new.power());
    assert_eq!(old.cells(), new.cells());
    let old = old
        .advance_to(&initial, 0.5, config(0.125), || Ok(()))
        .unwrap();
    let new = new
        .advance_to(&initial, 0.5, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(old.step.thermal.cells(), new.step.thermal.cells());
    assert_eq!(old.step.acoustic, new.step.acoustic);
    assert_eq!(old.bin_energy, new.bin_energy);
    assert_eq!(old.work, new.work);
    assert_eq!(old.heating.gas_states(), new.heating.gas_states());

    let uncorrelated = PorousCovarianceBin {
        covariance: covariance(9000.0, 4000.0, C64::ZERO),
        ..bins()[0]
    };
    let combined = PorousCrossSpectrumHeating::try_new(
        &initial,
        &routes(&materials),
        gas(),
        &[uncorrelated],
        options(),
        || Ok(()),
    )
    .unwrap();
    let individual: Vec<_> = (0..2)
        .map(|side| {
            let mut p = [C64::ZERO; 2];
            p[side] = C64::from_re((2.0_f64 * [9000.0, 4000.0][side]).sqrt());
            TemperaturePorousHeating::try_new(
                &initial,
                &routes(&materials),
                gas(),
                PorousHarmonicDrive {
                    angular_frequency_rad_s: uncorrelated.angular_frequency_rad_s,
                    area_m2: AREA,
                    reference_impedance_pa_s_m3: reference(),
                    incident_pressure_pa: p,
                },
                options().composition,
                || Ok(()),
            )
            .unwrap()
        })
        .collect();
    for i in 0..2 {
        close(
            combined.cells()[i].dissipation.viscous_power_w,
            individual
                .iter()
                .map(|h| h.cells()[i].dissipation.viscous_power_w)
                .sum(),
            2e-14,
        );
        close(
            combined.cells()[i].dissipation.thermal_power_w,
            individual
                .iter()
                .map(|h| h.cells()[i].dissipation.thermal_power_w)
                .sum(),
            2e-14,
        );
    }
    close(
        combined.power().incident_power_w,
        13000.0 / reference(),
        2e-15,
    );
    let zero = bins().map(|b| PorousCovarianceBin {
        covariance: PressureCovariance::try_new([C64::ZERO; 4]).unwrap(),
        ..b
    });
    let native = initial.advance_to(0.5, config(0.125), || Ok(())).unwrap();
    let source = PorousCrossSpectrumHeating::try_new(
        &initial,
        &routes(&materials),
        gas(),
        &zero,
        options(),
        || Ok(()),
    )
    .unwrap();
    let result = source
        .advance_to(&initial, 0.5, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(native.0.cells(), result.step.thermal.cells());
    assert_eq!(native.1, result.step.heat_report);
    assert_eq!(result.step.acoustic.deposited_j, 0.0);
}

#[test]
fn g3_partial_coherence_shares_gas_and_composition_and_preserves_events_latent_heat_and_continuation()
 {
    let materials = materials();
    let mut schedule = events();
    let Event::Birth(first) = &schedule[0].event else {
        panic!("birth")
    };
    schedule.push(TimedEvent {
        time_s: 0.25,
        event: Event::Boundary(
            IDS[0],
            Boundary {
                source_w: 0.4,
                ..first.boundary
            },
        ),
    });
    let initial = thermal(schedule, [&PLATEAU, &SENSIBLE]);
    let source = PorousCrossSpectrumHeating::try_new(
        &initial,
        &routes(&materials),
        gas(),
        &bins(),
        options(),
        || Ok(()),
    )
    .unwrap();
    assert!(source.bins()[0].angular_frequency_rad_s < source.bins()[1].angular_frequency_rad_s);
    let first = source
        .advance_to(&initial, 0.25, config(0.125), || Ok(()))
        .unwrap();
    let second = first
        .heating
        .advance_to(&first.step.thermal, 0.5, config(0.125), || Ok(()))
        .unwrap();
    let whole = source
        .advance_to(&initial, 0.5, config(0.125), || Ok(()))
        .unwrap();
    // Four accepted endpoints, one changed gas cell, TWO frequencies. The two
    // covariance components at each frequency do not double composition work.
    assert_eq!(whole.work.gas_evaluations, 4);
    assert_eq!(whole.work.stack_solves, 8);
    assert_eq!(second.step.thermal.cells(), whole.step.thermal.cells());
    assert_eq!(second.heating.power(), whole.heating.power());
    assert_eq!(whole.step.thermal.temperature_k(IDS[0]).unwrap(), 300.0);
    assert!(whole.step.thermal.cells()[&IDS[0]].specific_h_j_kg > 0.0);
    assert_eq!(whole.step.thermal.links(), initial.links());
    close(
        whole.step.heat_report.external_heat_j,
        whole.step.acoustic.deposited_j + 0.125,
        2e-13,
    );
    close(
        whole.step.acoustic.incident_j,
        0.5 * source.power().incident_power_w,
        2e-15,
    );
    close(
        whole.bin_energy.iter().map(|e| e.deposited_j).sum(),
        whole.step.acoustic.deposited_j,
        2e-15,
    );
    for j in 0..2 {
        close(
            first.bin_energy[j].outgoing_j + second.bin_energy[j].outgoing_j,
            whole.bin_energy[j].outgoing_j,
            2e-15,
        );
        close(
            first.bin_energy[j].deposited_j + second.bin_energy[j].deposited_j,
            whole.bin_energy[j].deposited_j,
            2e-15,
        );
        assert!(
            whole
                .heating
                .bin_power(j)
                .unwrap()
                .dissipation
                .absorbed_power_w
                > 0.0
        );
        assert!(
            whole
                .heating
                .bin_cells(j)
                .unwrap()
                .iter()
                .all(|c| c.dissipation.absorbed_power_w > 0.0)
        );
    }
    let plateau = thermal(events(), [&PLATEAU; 2]);
    let source = PorousCrossSpectrumHeating::try_new(
        &plateau,
        &routes(&materials),
        gas(),
        &bins(),
        options(),
        || Ok(()),
    )
    .unwrap();
    let accepted = source
        .advance_to(&plateau, 0.5, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(accepted.work, PorousFeedbackWork::default());
    close(
        accepted.step.acoustic.deposited_j,
        0.5 * source.power().dissipation.absorbed_power_w,
        3e-15,
    );
}

#[test]
fn g0_covariance_spectrum_refuses_ambiguous_frequencies_bad_bindings_and_late_candidates() {
    let materials = materials();
    let initial = thermal(events(), [&SENSIBLE; 2]);
    let routes = routes(&materials);
    let source = PorousCrossSpectrumHeating::try_new(
        &initial,
        &routes,
        gas(),
        &bins(),
        options(),
        || Ok(()),
    )
    .unwrap();
    let duplicate = [bins()[0], bins()[0]];
    assert!(matches!(
        PorousCrossSpectrumHeating::try_new(
            &initial,
            &routes,
            gas(),
            &duplicate,
            options(),
            || Ok(())
        ),
        Err(PorousHeatingError::Invalid(_))
    ));
    for mut bad in [bins(); 3].into_iter().enumerate() {
        match bad.0 {
            0 => bad.1[0].angular_frequency_rad_s = 0.0,
            1 => bad.1[0].angular_frequency_rad_s = core::f64::consts::TAU * 5000.0,
            _ => bad.1[1].area_m2 *= 2.0,
        }
        assert!(
            PorousCrossSpectrumHeating::try_new(
                &initial,
                &routes,
                gas(),
                &bad.1,
                options(),
                || Ok(())
            )
            .is_err()
        );
    }
    let small = PorousSpectrumOptions {
        maximum_cell_bins: 3,
        ..options()
    };
    assert!(matches!(
        PorousCrossSpectrumHeating::try_new(&initial, &routes, gas(), &bins(), small, || Ok(())),
        Err(PorousHeatingError::Thermal(DepositionError::Budget(_)))
    ));
    let narrow = IsobaricPorousGas::try_new(&SPEC, 101325.0, [299.0, 300.000001]).unwrap();
    let limited =
        PorousCrossSpectrumHeating::try_new(&initial, &routes, narrow, &bins(), options(), || {
            Ok(())
        })
        .unwrap();
    assert!(
        limited
            .advance_to(&initial, 0.125, config(0.125), || Ok(()))
            .is_err()
    );
    let mut changed = events();
    let Event::Birth(cell) = &mut changed[0].event else {
        panic!("birth")
    };
    cell.mass_kg *= 1.1;
    let stale = thermal(changed, [&SENSIBLE; 2]);
    assert!(
        source
            .advance_to(&stale, 0.125, config(0.125), || Ok(()))
            .is_err()
    );
    let law_change = thermal(events(), [&PLATEAU, &SENSIBLE]);
    assert!(
        source
            .advance_to(&law_change, 0.125, config(0.125), || Ok(()))
            .is_err()
    );
    let mut calls = 0;
    assert!(matches!(
        source.advance_to(&initial, 1.0, config(0.125), || {
            calls += 1;
            if calls > 180 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(PorousHeatingError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert!(initial.cells().values().all(|c| c.specific_h_j_kg == 0.0));
    assert!(source.gas_states().iter().all(|s| s.temperature == 300.0));
    let mut permuted = bins();
    permuted.reverse();
    let same =
        PorousCrossSpectrumHeating::try_new(&initial, &routes, gas(), &permuted, options(), || {
            Ok(())
        })
        .unwrap();
    assert_eq!(same.power(), source.power());
}
