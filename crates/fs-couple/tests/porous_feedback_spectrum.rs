//! Spectral power, native storage and coherent-field checks on authored
//! numerical fixtures; these are not measurements or material qualification.
use fs_couple::porous_heating::PorousHeatingError;
use fs_couple::porous_heating::feedback::spectrum::*;
use fs_couple::porous_heating::feedback::*;
use fs_duct::compose::ComposeOptions;
use fs_duct::porous::{PorousLayer, PorousLayerExcitation, PorousLayerStack};
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_math::c64::C64;
use fs_thermal::enthalpy::*;

#[path = "porous_feedback_spectrum/refinement.rs"]
mod refinement;

const AREA: f64 = 0.003_906_25;
const LENGTHS: [f64; 2] = [0.031_25, 0.046_875];
const CP: f64 = 200.0;
const IDS: [CellId; 2] = [CellId(19), CellId(3)];
static SPEC: GasSpec = GasSpec::dry_air_ussa1976();
struct Caloric(bool);
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if self.0 {
            if (0.0..=1e4).contains(&h) {
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
        Some([70 + u8::from(self.0); 32])
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
        maximum_bins: 128,
        maximum_cell_bins: 512,
        composition: ComposeOptions {
            maximum_blocks: 8,
            ..ComposeOptions::default()
        },
    }
}
fn materials() -> [PorousFluidLaw; 2] {
    std::array::from_fn(|i| {
        let mut card = PorousCard::melamine_uf_nguyen2024();
        card.name = "authored spectral absorption fixture";
        card.provenance = "numerical stationary-spectrum test, not measured qualification";
        card.flow_resistivity *= [1.0, 2.3][i];
        PorousFluidLaw::try_new(card, 1.81e-5).unwrap()
    })
}
fn routes(materials: &[PorousFluidLaw; 2]) -> [PorousTemperatureCell<'_>; 2] {
    std::array::from_fn(|i| PorousTemperatureCell {
        cell: IDS[i],
        material: &materials[i],
    })
}
fn gas() -> IsobaricPorousGas<'static> {
    IsobaricPorousGas::try_new(&SPEC, 101325.0, [280.0, 450.0]).unwrap()
}
fn drive(frequency_hz: f64, pressure: [C64; 2]) -> PorousHarmonicDrive {
    PorousHarmonicDrive {
        angular_frequency_rad_s: core::f64::consts::TAU * frequency_hz,
        area_m2: AREA,
        reference_impedance_pa_s_m3: GasState::try_new(&SPEC, 300.0, 101325.0)
            .unwrap()
            .characteristic_impedance
            / AREA,
        incident_pressure_pa: pressure,
    }
}
fn bins() -> [PorousHarmonicDrive; 3] {
    [
        drive(1900.0, [C64::new(90.0, -40.0), C64::new(-25.0, 50.0)]),
        drive(500.0, [C64::from_re(220.0), C64::new(70.0, -30.0)]),
        drive(1100.0, [C64::new(-50.0, 130.0), C64::from_re(80.0)]),
    ]
}
fn events(h: [f64; 2], sources: [f64; 2]) -> Vec<TimedEvent> {
    let mut start = 0.0;
    (0..2)
        .map(|i| {
            let end = start + LENGTHS[i];
            let event = TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: IDS[i],
                    segment_id: i as u64,
                    body_id: 72,
                    material_id: MaterialId(i as u64),
                    mass_kg: 100.0 * AREA * LENGTHS[i],
                    specific_h_j_kg: h[i],
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
                        source_w: sources[i],
                    },
                }),
            };
            start = end;
            event
        })
        .collect()
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
fn layers<'a>(
    materials: &'a [PorousFluidLaw; 2],
    temperatures: [f64; 2],
    drive: PorousHarmonicDrive,
) -> [PorousLayer<'a>; 2] {
    std::array::from_fn(|i| {
        let state = GasState::try_new(&SPEC, temperatures[i], 101325.0).unwrap();
        PorousLayer::try_new(
            &materials[i],
            &state,
            drive.angular_frequency_rad_s,
            LENGTHS[i],
            AREA,
            drive.reference_impedance_pa_s_m3,
        )
        .unwrap()
    })
}

// Independent p/v standing-wave reconstruction and Simpson volume integration,
// bypassing the production analytic integrals. Limited to these moderate slabs.
fn quadrature(field: &PorousLayerExcitation<'_, '_>) -> [f64; 2] {
    let response = field.layer().material_response();
    let (z, k) = response.characteristic();
    let (z, k) = (z.conj(), k.conj());
    let boundary = field.boundary_power();
    let p = boundary.pressure_pa[0];
    let v = boundary.inward_volume_flow_m3_s[0].scale(1.0 / AREA);
    let forward = (p + z * v).scale(0.5);
    let backward = (p - z * v).scale(0.5);
    let weights = [
        -0.5 * response.angular_frequency_rad_s() * response.dynamic_density_kg_m3().im,
        -0.5 * response.angular_frequency_rad_s() * response.dynamic_bulk_modulus_pa().recip().im,
    ];
    let dx = field.layer().thickness_m() / 512.0;
    let mut sum = [0.0; 2];
    for i in 0..=512 {
        let x = i as f64 * dx;
        let phase = C64::new((k.re * x).cos(), (k.re * x).sin()).scale((-k.im * x).exp());
        let (a, b) = (forward * phase, backward / phase);
        let magnitude = [((a - b) / z).norm_sq(), (a + b).norm_sq()];
        let w = if i == 0 || i == 512 {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for j in 0..2 {
            sum[j] += w * weights[j] * magnitude[j];
        }
    }
    sum.map(|value| value * AREA * dx / 3.0)
}

#[test]
fn g1_spectral_heat_matches_independent_volume_integrals_and_retains_bin_coherence() {
    let materials = materials();
    let initial = thermal(events([0.0, 25.0 * CP], [0.0; 2]), [&SENSIBLE; 2]);
    let drives = bins();
    let source = PorousSpectrumHeating::try_new(
        &initial,
        &routes(&materials),
        gas(),
        &drives,
        options(),
        || Ok(()),
    )
    .unwrap();
    let mut expected = [[0.0; 2]; 2];
    let mut exterior = [0.0; 2];
    for drive in drives {
        let slabs = layers(&materials, [300.0, 325.0], drive);
        let stack = PorousLayerStack::try_new(&slabs, options().composition, || false).unwrap();
        let field = stack
            .prepare_excitation(drive.incident_pressure_pa, || false)
            .unwrap();
        for i in 0..2 {
            let work = quadrature(&field.layer_excitations()[i]);
            for j in 0..2 {
                expected[i][j] += work[j];
            }
        }
        exterior[0] += field.boundary_power().incident_power_w;
        exterior[1] += field.boundary_power().outgoing_power_w;
    }
    for i in 0..2 {
        close(
            source.cells()[i].dissipation.viscous_power_w,
            expected[i][0],
            3e-9,
        );
        close(
            source.cells()[i].dissipation.thermal_power_w,
            expected[i][1],
            3e-9,
        );
    }
    close(source.power().incident_power_w, exterior[0], 2e-15);
    close(source.power().outgoing_power_w, exterior[1], 2e-15);
    close(source.power().residual_w, 0.0, 0.0);
    assert!(
        source
            .bins()
            .windows(2)
            .all(|p| p[0].drive().angular_frequency_rad_s < p[1].drive().angular_frequency_rad_s)
    );
    let mut permuted = drives;
    permuted.swap(0, 2);
    let reordered = PorousSpectrumHeating::try_new(
        &initial,
        &routes(&materials),
        gas(),
        &permuted,
        options(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(source.cells(), reordered.cells());
    assert_eq!(source.power(), reordered.power());

    // Overall phase of one independent bin has no observable thermal effect.
    let mut rotated = drives;
    rotated[1].incident_pressure_pa = rotated[1]
        .incident_pressure_pa
        .map(|p| C64::new(-p.im, p.re));
    let rotated = PorousSpectrumHeating::try_new(
        &initial,
        &routes(&materials),
        gas(),
        &rotated,
        options(),
        || Ok(()),
    )
    .unwrap();
    close(
        source.power().dissipation.absorbed_power_w,
        rotated.power().dissipation.absorbed_power_w,
        2e-14,
    );
    // Relative phase WITHIN one bin changes its coherent standing-wave work.
    let mut opposite = drives;
    opposite[1].incident_pressure_pa[1] = -opposite[1].incident_pressure_pa[1];
    let opposite = PorousSpectrumHeating::try_new(
        &initial,
        &routes(&materials),
        gas(),
        &opposite,
        options(),
        || Ok(()),
    )
    .unwrap();
    assert!(
        (source.power().dissipation.absorbed_power_w
            - opposite.power().dissipation.absorbed_power_w)
            .abs()
            > 1e-4
    );
}

#[test]
fn g1_one_sided_angular_density_has_correct_power_and_frequency_quadrature_refines() {
    let materials = materials();
    let initial = thermal(events([0.0; 2], [0.0; 2]), [&PLATEAU; 2]);
    let density = [C64::new(0.5, 0.2), C64::new(-0.1, 0.3)];
    let template = drive(1000.0, [C64::ZERO; 2]);
    let r = template.reference_impedance_pa_s_m3;
    let omega0 = core::f64::consts::TAU * 200.0;
    let width = core::f64::consts::TAU * 2400.0;
    let bin =
        angular_pressure_density_bin(template.angular_frequency_rad_s, AREA, r, density, 123.0)
            .unwrap();
    let source = PorousSpectrumHeating::try_new(
        &initial,
        &routes(&materials),
        gas(),
        &[bin],
        options(),
        || Ok(()),
    )
    .unwrap();
    close(
        source.power().incident_power_w,
        123.0 * density.iter().map(|p| p.norm_sq()).sum::<f64>() / r,
        2e-15,
    );
    let from_hz = angular_pressure_density_bin(
        template.angular_frequency_rad_s,
        AREA,
        r,
        density.map(|p| p.scale(1.0 / core::f64::consts::TAU.sqrt())),
        123.0 * core::f64::consts::TAU,
    )
    .unwrap();
    for j in 0..2 {
        close(
            bin.incident_pressure_pa[j].re,
            from_hz.incident_pressure_pa[j].re,
            2e-15,
        );
    }

    // Direct original constitutive/stack integration in omega is independent
    // of the adapter's bin conversion and spectral accumulation.
    let mut reference = 0.0;
    for i in 0..=1024 {
        let omega = omega0 + width * i as f64 / 1024.0;
        let drive = PorousHarmonicDrive {
            angular_frequency_rad_s: omega,
            incident_pressure_pa: density.map(|p| p.scale(core::f64::consts::SQRT_2)),
            ..template
        };
        let slabs = layers(&materials, [300.0; 2], drive);
        let stack = PorousLayerStack::try_new(&slabs, options().composition, || false).unwrap();
        let work = stack
            .prepare_excitation(drive.incident_pressure_pa, || false)
            .unwrap()
            .dissipation()
            .absorbed_power_w;
        let w = if i == 0 || i == 1024 {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        reference += w * work;
    }
    reference *= width / 1024.0 / 3.0;
    let mut previous = f64::INFINITY;
    for count in [8, 16, 32] {
        let dw = width / count as f64;
        let drives: Vec<_> = (0..count)
            .map(|i| {
                angular_pressure_density_bin(omega0 + (i as f64 + 0.5) * dw, AREA, r, density, dw)
                    .unwrap()
            })
            .collect();
        let source = PorousSpectrumHeating::try_new(
            &initial,
            &routes(&materials),
            gas(),
            &drives,
            options(),
            || Ok(()),
        )
        .unwrap();
        close(
            source.power().incident_power_w,
            width * density.iter().map(|p| p.norm_sq()).sum::<f64>() / r,
            3e-15,
        );
        let error = (source.power().dissipation.absorbed_power_w - reference).abs();
        assert!(
            error < 0.4 * previous,
            "spectral quadrature error {error} after {previous}"
        );
        previous = error;
    }
    assert!(previous < 5e-4 * reference);
}

#[test]
fn g1_shared_gas_cache_latent_storage_zero_drive_and_single_bin_parity() {
    let materials = materials();
    let initial = thermal(events([0.0; 2], [0.1; 2]), [&PLATEAU; 2]);
    let source = PorousSpectrumHeating::try_new(
        &initial,
        &routes(&materials),
        gas(),
        &bins(),
        options(),
        || Ok(()),
    )
    .unwrap();
    let plateau = source
        .advance_to(&initial, 1.0, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(plateau.work, PorousFeedbackWork::default());
    for i in 0..2 {
        close(
            plateau.step.cell_energy_j[i].1,
            source.cells()[i].dissipation.absorbed_power_w,
            2e-15,
        );
        let cell = &plateau.step.thermal.cells()[&IDS[i]];
        close(
            cell.mass_kg * cell.specific_h_j_kg,
            plateau.step.cell_energy_j[i].1 + 0.1,
            2e-14,
        );
        assert_eq!(plateau.step.thermal.temperature_k(IDS[i]).unwrap(), 300.0);
    }
    close(
        plateau.bin_energy.iter().map(|e| e.deposited_j).sum(),
        plateau.step.acoustic.deposited_j,
        3e-15,
    );
    close(
        plateau.bin_energy.iter().map(|e| e.outgoing_j).sum(),
        plateau.step.acoustic.outgoing_j,
        3e-15,
    );

    let selective = thermal(events([0.0; 2], [0.0; 2]), [&PLATEAU, &SENSIBLE]);
    let shared = PorousSpectrumHeating::try_new(
        &selective,
        &routes(&materials),
        gas(),
        &bins(),
        options(),
        || Ok(()),
    )
    .unwrap();
    let next = shared
        .advance_to(&selective, 0.125, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(next.work.gas_evaluations, 1);
    assert_eq!(next.work.stack_solves, 3);
    close(
        next.step.acoustic.deposited_j,
        0.125 * shared.power().dissipation.absorbed_power_w,
        2e-15,
    );
    for bin in next.heating.bins() {
        assert_eq!(bin.gas_states(), next.heating.gas_states());
    }

    let initial = thermal(events([0.0, 10.0 * CP], [0.2, -0.1]), [&SENSIBLE; 2]);
    let single = TemperaturePorousHeating::try_new(
        &initial,
        &routes(&materials),
        gas(),
        bins()[0],
        options().composition,
        || Ok(()),
    )
    .unwrap();
    let spectrum = PorousSpectrumHeating::try_new(
        &initial,
        &routes(&materials),
        gas(),
        &[bins()[0]],
        options(),
        || Ok(()),
    )
    .unwrap();
    let one = single
        .advance_to(&initial, 1.0, config(0.125), || Ok(()))
        .unwrap();
    let many = spectrum
        .advance_to(&initial, 1.0, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(one.step.thermal.cells(), many.step.thermal.cells());
    assert_eq!(one.step.acoustic, many.step.acoustic);
    assert_eq!(many.bin_energy[0], many.step.acoustic);
    assert_eq!(one.work, many.work);
    assert_eq!(one.heating.gas_states(), many.heating.gas_states());

    let zero_bins = bins().map(|b| PorousHarmonicDrive {
        incident_pressure_pa: [C64::ZERO; 2],
        ..b
    });
    let zero = PorousSpectrumHeating::try_new(
        &initial,
        &routes(&materials),
        gas(),
        &zero_bins,
        options(),
        || Ok(()),
    )
    .unwrap();
    let native = initial.advance_to(1.0, config(0.125), || Ok(())).unwrap();
    let coupled = zero
        .advance_to(&initial, 1.0, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(native.0.cells(), coupled.step.thermal.cells());
    assert_eq!(native.1, coupled.step.heat_report);
    assert_eq!(coupled.step.acoustic.deposited_j, 0.0);
}

#[test]
fn g3_spectral_events_continuation_domains_and_late_cancellation_are_atomic() {
    let materials = materials();
    let mut schedule = events([0.0; 2], [0.1, -0.05]);
    schedule.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 1,
            a: IDS[0],
            b: IDS[1],
            area_m2: AREA,
            conductance_w_k: 0.02,
        }),
    });
    let Event::Birth(first) = &schedule[0].event else {
        panic!("birth")
    };
    schedule.push(TimedEvent {
        time_s: 0.25,
        event: Event::Boundary(
            IDS[0],
            Boundary {
                source_w: 0.3,
                ..first.boundary
            },
        ),
    });
    let initial = thermal(schedule.clone(), [&SENSIBLE; 2]);
    let routes = routes(&materials);
    let source =
        PorousSpectrumHeating::try_new(&initial, &routes, gas(), &bins(), options(), || Ok(()))
            .unwrap();
    let first = source
        .advance_to(&initial, 0.5, config(0.125), || Ok(()))
        .unwrap();
    let second = first
        .heating
        .advance_to(&first.step.thermal, 1.0, config(0.125), || Ok(()))
        .unwrap();
    let whole = source
        .advance_to(&initial, 1.0, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(second.step.thermal.cells(), whole.step.thermal.cells());
    assert_eq!(second.heating.power(), whole.heating.power());
    for j in 0..3 {
        close(
            first.bin_energy[j].deposited_j + second.bin_energy[j].deposited_j,
            whole.bin_energy[j].deposited_j,
            3e-15,
        );
        close(
            first.bin_energy[j].outgoing_j + second.bin_energy[j].outgoing_j,
            whole.bin_energy[j].outgoing_j,
            3e-15,
        );
    }
    close(
        whole.step.heat_report.external_heat_j,
        whole.step.acoustic.deposited_j + 0.05 * 0.25 + 0.25 * 0.75,
        3e-14,
    );
    close(
        whole.step.acoustic.incident_j,
        source.power().incident_power_w,
        3e-15,
    );
    assert_eq!(whole.step.thermal.links(), initial.links());
    assert_eq!(whole.step.thermal.cells()[&IDS[0]].boundary.source_w, 0.3);

    for limited in [
        PorousSpectrumOptions {
            maximum_bins: 2,
            ..options()
        },
        PorousSpectrumOptions {
            maximum_cell_bins: 5,
            ..options()
        },
    ] {
        assert!(matches!(
            PorousSpectrumHeating::try_new(&initial, &routes, gas(), &bins(), limited, || Ok(())),
            Err(PorousHeatingError::Thermal(DepositionError::Budget(_)))
        ));
    }
    assert!(
        PorousSpectrumHeating::try_new(
            &initial,
            &routes,
            gas(),
            &[bins()[0], bins()[0]],
            options(),
            || Ok(())
        )
        .is_err()
    );
    for invalid in [0.0, -1.0, f64::NAN, core::f64::consts::TAU * 5000.0] {
        let mut bad = bins();
        bad[2].angular_frequency_rad_s = invalid;
        assert!(
            PorousSpectrumHeating::try_new(&initial, &routes, gas(), &bad, options(), || Ok(()))
                .is_err()
        );
    }
    let mut wrong_area = bins();
    wrong_area[1].area_m2 *= 2.0;
    assert!(
        PorousSpectrumHeating::try_new(&initial, &routes, gas(), &wrong_area, options(), || Ok(()))
            .is_err()
    );
    for weight in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(angular_pressure_density_bin(1000.0, AREA, 1e5, [C64::ONE; 2], weight).is_err());
    }
    let narrow = IsobaricPorousGas::try_new(&SPEC, 101325.0, [299.0, 300.000001]).unwrap();
    let limited =
        PorousSpectrumHeating::try_new(&initial, &routes, narrow, &bins(), options(), || Ok(()))
            .unwrap();
    assert!(
        limited
            .advance_to(&initial, 0.125, config(0.125), || Ok(()))
            .is_err()
    );
    let Event::Birth(cell) = &mut schedule[0].event else {
        panic!("birth")
    };
    cell.mass_kg *= 1.1;
    let stale = thermal(schedule, [&SENSIBLE; 2]);
    assert!(
        source
            .advance_to(&stale, 0.125, config(0.125), || Ok(()))
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        source.advance_to(&initial, 1.0, config(0.125), || {
            polls += 1;
            if polls > 180 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(PorousHeatingError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert!(
        initial
            .cells()
            .values()
            .all(|cell| cell.specific_h_j_kg == 0.0)
    );
    assert!(
        source
            .gas_states()
            .iter()
            .all(|state| state.temperature == 300.0)
    );
    assert_eq!(
        source.power(),
        PorousSpectrumHeating::try_new(&initial, &routes, gas(), &bins(), options(), || Ok(()))
            .unwrap()
            .power()
    );
}
