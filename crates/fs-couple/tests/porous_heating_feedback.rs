//! Native gas/porous feedback on authored lumped caloric fixtures. These checks
//! establish numerical/model composition, not qualification of a real specimen.
use fs_couple::porous_heating::feedback::*;
use fs_couple::porous_heating::{PorousHeatCell, PorousHeatingError, PorousStackHeating};
use fs_duct::compose::ComposeOptions;
use fs_duct::porous::{PorousLayer, PorousLayerStack};
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_math::c64::C64;
use fs_thermal::enthalpy::*;

#[path = "porous_heating_feedback/refinement.rs"]
mod refinement;
#[path = "porous_heating_feedback/termination.rs"]
mod termination;

const AREA: f64 = 0.004;
const CP: f64 = 200.0;
static SPEC: GasSpec = GasSpec::dry_air_ussa1976();
struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h / CP;
        if t.is_finite() && t > 0.0 {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / CP
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([40; 32])
    }
}
struct Plateau;
impl ThermalLaw for Plateau {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if (0.0..=1e4).contains(&h) {
            Ok(300.0)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.0
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([41; 32])
    }
}
static SENSIBLE: Sensible = Sensible;
static PLATEAU: Plateau = Plateau;
fn id(index: usize) -> CellId {
    CellId(30 - 7 * index as u64)
}
fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-9,
    }
}
fn options() -> ComposeOptions {
    ComposeOptions {
        maximum_blocks: 16,
        ..ComposeOptions::default()
    }
}
fn materials() -> [PorousFluidLaw; 2] {
    std::array::from_fn(|i| {
        let mut card = PorousCard::melamine_uf_nguyen2024();
        card.name = "authored porous feedback fixture";
        card.provenance = "numerical thermal feedback law, not measured qualification";
        card.flow_resistivity *= [1.0, 2.3][i];
        PorousFluidLaw::try_new(card, 1.81e-5).unwrap()
    })
}
fn cells(materials: &[PorousFluidLaw; 2]) -> [PorousTemperatureCell<'_>; 2] {
    std::array::from_fn(|i| PorousTemperatureCell {
        cell: id(i),
        material: &materials[i],
    })
}
fn gas() -> IsobaricPorousGas<'static> {
    IsobaricPorousGas::try_new(&SPEC, 101325.0, [280.0, 450.0]).unwrap()
}
fn drive() -> PorousHarmonicDrive {
    PorousHarmonicDrive {
        angular_frequency_rad_s: core::f64::consts::TAU * 1100.0,
        area_m2: AREA,
        reference_impedance_pa_s_m3: GasState::try_new(&SPEC, 300.0, 101325.0)
            .unwrap()
            .characteristic_impedance
            / AREA,
        incident_pressure_pa: [C64::from_re(220.0), C64::new(70.0, -30.0)],
    }
}
fn events(h: [f64; 2], sources: [f64; 2]) -> Vec<TimedEvent> {
    let mut result = Vec::new();
    let ends = [0.0, 0.03, 0.08];
    for i in 0..2 {
        let length = ends[i + 1] - ends[i];
        result.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: id(i),
                segment_id: i as u64,
                body_id: 90,
                material_id: MaterialId(i as u64),
                mass_kg: 100.0 * AREA * length,
                specific_h_j_kg: h[i],
                geometry: CellGeometry {
                    start_m: [ends[i], 0.0, 0.0],
                    end_m: [ends[i + 1], 0.0, 0.0],
                    volume_m3: AREA * length,
                    surface_area_m2: 2.0 * AREA + 4.0 * AREA.sqrt() * length,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: sources[i],
                },
            }),
        });
    }
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
fn native_layers<'a>(
    materials: &'a [PorousFluidLaw; 2],
    temperatures: [f64; 2],
) -> Vec<PorousLayer<'a>> {
    let drive = drive();
    let ends = [0.0, 0.03, 0.08];
    (0..2)
        .map(|i| {
            let state = GasState::try_new(&SPEC, temperatures[i], 101325.0).unwrap();
            PorousLayer::try_new(
                &materials[i],
                &state,
                drive.angular_frequency_rad_s,
                ends[i + 1] - ends[i],
                AREA,
                drive.reference_impedance_pa_s_m3,
            )
            .unwrap()
        })
        .collect()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

#[test]
fn g1_cell_temperature_gradient_uses_original_isobaric_gas_and_coherent_loss() {
    let materials = materials();
    let thermal = thermal(events([0.0, CP * 40.0], [0.0; 2]), [&SENSIBLE; 2]);
    let heating = TemperaturePorousHeating::try_new(
        &thermal,
        &cells(&materials),
        gas(),
        drive(),
        options(),
        || Ok(()),
    )
    .unwrap();
    let original_layers = native_layers(&materials, [300.0, 340.0]);
    let original_stack = PorousLayerStack::try_new(&original_layers, options(), || false).unwrap();
    let original = original_stack
        .prepare_excitation(drive().incident_pressure_pa, || false)
        .unwrap();
    for i in 0..2 {
        let expected =
            GasState::try_new(&SPEC, thermal.temperature_k(id(i)).unwrap(), 101325.0).unwrap();
        assert_eq!(heating.gas_states()[i], expected);
        assert_eq!(
            *heating.acoustic_layers()[i].material_response(),
            *original_layers[i].material_response()
        );
        assert!(std::ptr::eq(
            heating.acoustic_layers()[i].material(),
            &materials[i]
        ));
        let loss = original.layer_excitations()[i]
            .dissipation_between(0.0, original_layers[i].thickness_m())
            .unwrap();
        assert_eq!(heating.cells()[i].dissipation, loss);
    }
    // Independent perfect-gas identities: density falls and c rises together;
    // the original viscous resistivity also follows its own gas viscosity.
    let states = heating.gas_states();
    close(states[1].density / states[0].density, 300.0 / 340.0, 2e-15);
    close(
        states[1].sound_speed / states[0].sound_speed,
        (340.0_f64 / 300.0).sqrt(),
        2e-15,
    );
    assert!(states[1].dynamic_viscosity > states[0].dynamic_viscosity);
    assert_eq!(heating.dissipation(), original.dissipation());
    let uniform_layers = native_layers(&materials, [320.0; 2]);
    let uniform_stack = PorousLayerStack::try_new(&uniform_layers, options(), || false).unwrap();
    let uniform = uniform_stack
        .prepare_excitation(drive().incident_pressure_pa, || false)
        .unwrap();
    assert!(
        (heating.dissipation().absorbed_power_w - uniform.dissipation().absorbed_power_w).abs()
            > 1e-4 * heating.dissipation().absorbed_power_w
    );
    assert_eq!(
        heating.boundary_power().incident_power_w,
        uniform.boundary_power().incident_power_w
    );
}

#[test]
fn g1_latent_plateau_retains_native_heat_and_only_changed_temperatures_rebuild() {
    let materials = materials();
    let initial = thermal(events([0.0; 2], [0.03, 0.07]), [&PLATEAU; 2]);
    let heating = TemperaturePorousHeating::try_new(
        &initial,
        &cells(&materials),
        gas(),
        drive(),
        options(),
        || Ok(()),
    )
    .unwrap();
    let layers = native_layers(&materials, [300.0; 2]);
    let stack = PorousLayerStack::try_new(&layers, options(), || false).unwrap();
    let excitation = stack
        .prepare_excitation(drive().incident_pressure_pa, || false)
        .unwrap();
    let routes = [
        vec![PorousHeatCell {
            cell: id(0),
            end_depth_m: layers[0].thickness_m(),
        }],
        vec![PorousHeatCell {
            cell: id(1),
            end_depth_m: layers[1].thickness_m(),
        }],
    ];
    let frozen =
        PorousStackHeating::try_new(&excitation, &initial, &[&routes[0], &routes[1]], 2, || {
            Ok(())
        })
        .unwrap();
    let held = frozen
        .advance_to(&initial, 0.5, config(0.125), || Ok(()))
        .unwrap();
    let evolving = heating
        .advance_to(&initial, 0.5, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(evolving.step.thermal.cells(), held.thermal.cells());
    assert_eq!(evolving.step.acoustic, held.acoustic);
    assert_eq!(evolving.step.heat_report, held.heat_report);
    assert_eq!(evolving.work, PorousFeedbackWork::default());
    assert!(
        evolving
            .step
            .thermal
            .cells()
            .values()
            .all(|cell| cell.specific_h_j_kg > 0.0)
    );
    assert!(
        evolving
            .heating
            .gas_states()
            .iter()
            .all(|state| state.temperature == 300.0)
    );

    let selective = thermal(events([0.0; 2], [0.0; 2]), [&PLATEAU, &SENSIBLE]);
    let heating = TemperaturePorousHeating::try_new(
        &selective,
        &cells(&materials),
        gas(),
        drive(),
        options(),
        || Ok(()),
    )
    .unwrap();
    let before = heating.dissipation().absorbed_power_w;
    let next = heating
        .advance_to(&selective, 0.125, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(next.work.gas_evaluations, 1);
    assert_eq!(next.work.stack_solves, 1);
    assert_eq!(next.heating.gas_states()[0], heating.gas_states()[0]);
    assert!(next.heating.gas_states()[1].temperature > heating.gas_states()[1].temperature);
    close(next.step.acoustic.deposited_j, 0.125 * before, 2e-15);
    for i in 0..2 {
        assert_eq!(
            next.heating.gas_states()[i].temperature,
            next.step.thermal.temperature_k(id(i)).unwrap()
        );
    }
}

#[test]
fn g3_feedback_events_and_continuation_keep_fixed_exterior_power_and_one_credit() {
    let materials = materials();
    let mut schedule = events([0.0, 10.0 * CP], [0.1, -0.05]);
    schedule.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 1,
            a: id(0),
            b: id(1),
            area_m2: AREA,
            conductance_w_k: 0.01,
        }),
    });
    let Event::Birth(first_cell) = &schedule[0].event else {
        panic!("birth")
    };
    schedule.push(TimedEvent {
        time_s: 0.25,
        event: Event::Boundary(
            id(0),
            Boundary {
                source_w: 0.3,
                ..first_cell.boundary
            },
        ),
    });
    let initial = thermal(schedule, [&SENSIBLE; 2]);
    let heating = TemperaturePorousHeating::try_new(
        &initial,
        &cells(&materials),
        gas(),
        drive(),
        options(),
        || Ok(()),
    )
    .unwrap();
    let first = heating
        .advance_to(&initial, 0.5, config(0.125), || Ok(()))
        .unwrap();
    let second = first
        .heating
        .advance_to(&first.step.thermal, 1.0, config(0.125), || Ok(()))
        .unwrap();
    let whole = heating
        .advance_to(&initial, 1.0, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(second.step.thermal.cells(), whole.step.thermal.cells());
    assert_eq!(second.heating.gas_states(), whole.heating.gas_states());
    close(
        first.step.acoustic.deposited_j + second.step.acoustic.deposited_j,
        whole.step.acoustic.deposited_j,
        3e-14,
    );
    close(
        whole.step.heat_report.external_heat_j,
        whole.step.acoustic.deposited_j + 0.05 * 0.25 + 0.25 * 0.75,
        3e-13,
    );
    close(
        whole.step.acoustic.incident_j,
        heating.boundary_power().incident_power_w,
        3e-14,
    );
    close(whole.step.acoustic.residual_j, 0.0, 3e-13);
    assert_eq!(whole.step.thermal.links(), initial.links());
    assert_eq!(whole.step.thermal.cells()[&id(0)].boundary.source_w, 0.3);
    assert_eq!(initial.cells()[&id(0)].boundary.source_w, 0.1);
    assert_eq!(heating.gas_states()[0].temperature, 300.0);
}

#[test]
fn g0_temperature_domains_geometry_material_bindings_and_cancellation_remain_atomic() {
    let materials = materials();
    let schedule = events([0.0; 2], [0.0; 2]);
    let initial = thermal(schedule.clone(), [&SENSIBLE; 2]);
    let routes = cells(&materials);
    assert!(matches!(
        IsobaricPorousGas::try_new(&SPEC, 101325.0, [1.0, 300.0]),
        Err(PorousHeatingError::Material(_))
    ));
    assert!(matches!(
        IsobaricPorousGas::try_new(&SPEC, 1e8, [300.0, 310.0]),
        Err(PorousHeatingError::Material(_))
    ));
    assert!(IsobaricPorousGas::try_new(&SPEC, 101325.0, [310.0, 300.0]).is_err());
    let heating =
        TemperaturePorousHeating::try_new(&initial, &routes, gas(), drive(), options(), || Ok(()))
            .unwrap();
    let narrow = IsobaricPorousGas::try_new(&SPEC, 101325.0, [299.0, 300.000001]).unwrap();
    let limited =
        TemperaturePorousHeating::try_new(&initial, &routes, narrow, drive(), options(), || Ok(()))
            .unwrap();
    assert!(
        limited
            .advance_to(&initial, 0.125, config(0.125), || Ok(()))
            .is_err()
    );
    assert!(
        TemperaturePorousHeating::try_new(
            &initial,
            &routes,
            gas(),
            PorousHarmonicDrive {
                angular_frequency_rad_s: core::f64::consts::TAU * 5000.0,
                ..drive()
            },
            options(),
            || Ok(())
        )
        .is_err()
    );
    assert!(
        TemperaturePorousHeating::try_new(
            &initial,
            &[routes[0], routes[0]],
            gas(),
            drive(),
            options(),
            || Ok(())
        )
        .is_err()
    );
    assert!(
        TemperaturePorousHeating::try_new(
            &initial,
            &routes,
            gas(),
            drive(),
            ComposeOptions {
                maximum_blocks: 1,
                ..options()
            },
            || Ok(())
        )
        .is_err()
    );
    for volume in [true, false] {
        let mut changed = schedule.clone();
        let Event::Birth(cell) = &mut changed[0].event else {
            panic!("birth")
        };
        if volume {
            cell.geometry.volume_m3 *= 1.1;
        } else {
            cell.mass_kg *= 1.1;
        }
        let changed = thermal(changed, [&SENSIBLE; 2]);
        if volume {
            assert!(
                TemperaturePorousHeating::try_new(
                    &changed,
                    &routes,
                    gas(),
                    drive(),
                    options(),
                    || Ok(())
                )
                .is_err()
            );
        }
        assert!(
            heating
                .advance_to(&changed, 0.125, config(0.125), || Ok(()))
                .is_err()
        );
    }
    let changed_law = thermal(schedule, [&PLATEAU, &SENSIBLE]);
    assert!(
        heating
            .advance_to(&changed_law, 0.125, config(0.125), || Ok(()))
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        heating.advance_to(&initial, 1.0, config(0.125), || {
            polls += 1;
            if polls > 60 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(PorousHeatingError::Thermal(DepositionError::Cancelled))
    ));
    assert!(matches!(
        TemperaturePorousHeating::try_new(&initial, &routes, gas(), drive(), options(), || Err(
            DepositionError::Cancelled
        )),
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
        heating
            .gas_states()
            .iter()
            .all(|state| state.temperature == 300.0)
    );
}
