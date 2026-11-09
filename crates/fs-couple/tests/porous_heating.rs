//! Independent acoustic-volume and native finite-storage checks. Authored
//! numerical material fixtures establish implementation, not calibration.
use fs_couple::porous_heating::*;
use fs_duct::porous::{PorousLayer, PorousLayerExcitation};
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_math::c64::C64;
use fs_thermal::enthalpy::*;

#[path = "porous_heating/stack.rs"]
mod stack;
#[path = "porous_heating/transient.rs"]
mod transient;

const LENGTH: f64 = 0.08;
const AREA: f64 = 0.004;
const DENSITY: f64 = 100.0;
const CP: f64 = 200.0;

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
        Some([11; 32])
    }
}
struct Latent;
impl ThermalLaw for Latent {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=300.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if h < 100.0 {
            300.0 + h / CP
        } else if h <= 200.0 {
            300.5
        } else {
            300.5 + (h - 200.0) / CP
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / CP
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([22; 32])
    }
}
static SENSIBLE: Sensible = Sensible;
static LATENT: Latent = Latent;
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn material() -> PorousFluidLaw {
    PorousFluidLaw::try_new(PorousCard::melamine_uf_nguyen2024(), 1.81e-5).unwrap()
}
fn layer(material: &PorousFluidLaw) -> PorousLayer<'_> {
    let gas = gas();
    PorousLayer::try_new(
        material,
        &gas,
        core::f64::consts::TAU * 1100.0,
        LENGTH,
        AREA,
        gas.characteristic_impedance / AREA,
    )
    .unwrap()
}
fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-10,
    }
}
fn id(i: usize) -> CellId {
    CellId(1000 - 7 * i as u64)
}
fn point(depth: f64) -> [f64; 3] {
    [0.2 + 0.6 * depth, -0.3 + 0.8 * depth, 0.1]
}
fn routes(ends: &[f64]) -> Vec<PorousHeatCell> {
    ends.iter()
        .enumerate()
        .map(|(i, &end_depth_m)| PorousHeatCell {
            cell: id(i),
            end_depth_m,
        })
        .collect()
}
fn events(ends: &[f64], initial_h: &[f64], sources: &[f64], conductance: f64) -> Vec<TimedEvent> {
    let mut result = Vec::new();
    let mut start = 0.0;
    for (i, &end) in ends.iter().enumerate() {
        let width = end - start;
        result.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: id(i),
                segment_id: i as u64,
                body_id: 70,
                material_id: MaterialId(5),
                mass_kg: DENSITY * AREA * width,
                specific_h_j_kg: initial_h[i],
                geometry: CellGeometry {
                    start_m: point(start),
                    end_m: point(end),
                    volume_m3: AREA * width,
                    surface_area_m2: 2.0 * AREA + 4.0 * AREA.sqrt() * width,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: sources[i],
                },
            }),
        });
        start = end;
    }
    if conductance > 0.0 {
        for i in 1..ends.len() {
            result.push(TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: i as u64,
                    a: id(i - 1),
                    b: id(i),
                    area_m2: AREA,
                    conductance_w_k: conductance,
                }),
            });
        }
    }
    result
}
fn initial(events: Vec<TimedEvent>, law: &'static dyn ThermalLaw) -> DepositionSimulation<'static> {
    DepositionSimulation::new(
        events,
        &[(MaterialId(5), law)],
        DepositionBudget {
            max_cells: 32,
            max_links: 32,
            max_events: 64,
            max_substeps: 10000,
        },
    )
    .unwrap()
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

// Independent standing-wave reconstruction and Simpson volume integral. Unlike
// production, the backward wave is anchored at x=0 and may grow: the oracle is
// deliberately restricted to this moderately thick fixture.
fn quadrature(excitation: &PorousLayerExcitation<'_, '_>, start: f64, end: f64) -> [f64; 2] {
    let response = excitation.layer().material_response();
    let (z, k) = response.characteristic();
    let (z, k) = (z.conj(), k.conj());
    let boundary = excitation.boundary_power();
    let p = boundary.pressure_pa[0];
    let v = boundary.inward_volume_flow_m3_s[0].scale(1.0 / AREA);
    let forward = (p + z * v).scale(0.5);
    let backward = (p - z * v).scale(0.5);
    let weights = [
        -0.5 * response.angular_frequency_rad_s() * response.dynamic_density_kg_m3().im,
        -0.5 * response.angular_frequency_rad_s() * response.dynamic_bulk_modulus_pa().recip().im,
    ];
    let panels = 1024;
    let dx = (end - start) / panels as f64;
    let mut sum = [0.0; 2];
    for i in 0..=panels {
        let x = start + i as f64 * dx;
        let phase = C64::new((k.re * x).cos(), (k.re * x).sin()).scale((-k.im * x).exp());
        let a = forward * phase;
        let b = backward / phase;
        let fields = [((a - b) / z).norm_sq(), (a + b).norm_sq()];
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for j in 0..2 {
            sum[j] += weight * weights[j] * fields[j];
        }
    }
    sum.map(|value| value * AREA * dx / 3.0)
}

#[test]
fn g1_coherent_cell_heat_matches_volume_integral_and_partition_independent_energy() {
    let material = material();
    let layer = layer(&material);
    let excitation = layer
        .prepare_excitation([C64::from_re(100.0), C64::new(35.0, -20.0)])
        .unwrap();
    let mut totals = Vec::new();
    for n in [1, 3, 11] {
        let mut ends: Vec<_> = (1..=n).map(|i| LENGTH * i as f64 / n as f64).collect();
        *ends.last_mut().unwrap() = LENGTH;
        let thermal = initial(events(&ends, &vec![0.0; n], &vec![0.0; n], 0.0), &SENSIBLE);
        let heating =
            PorousHeating::try_new(&excitation, &thermal, &routes(&ends), 32, || Ok(())).unwrap();
        for cell in heating.cells() {
            let expected = quadrature(&excitation, cell.start_depth_m, cell.end_depth_m);
            close(
                cell.dissipation.viscous_power_w,
                expected[0],
                3e-11 * expected[0],
            );
            close(
                cell.dissipation.thermal_power_w,
                expected[1],
                3e-11 * expected[1],
            );
        }
        let step = heating
            .advance_to(&thermal, 1.3, config(0.17), || Ok(()))
            .unwrap();
        for ((cell, &(energy_id, energy)), &route) in heating
            .cells()
            .iter()
            .zip(&step.cell_energy_j)
            .zip(&routes(&ends))
        {
            assert_eq!(cell.cell, route.cell);
            assert_eq!(energy_id, route.cell);
            let original = &thermal.cells()[&route.cell];
            let expected = 1.3 * cell.dissipation.absorbed_power_w;
            close(energy, expected, 1e-14);
            close(
                step.thermal.cells()[&route.cell].specific_h_j_kg * original.mass_kg,
                expected,
                1e-14,
            );
            close(
                step.thermal.temperature_k(route.cell).unwrap(),
                300.0 + expected / (CP * original.mass_kg),
                1e-12,
            );
            assert_eq!(step.thermal.cells()[&route.cell].boundary.source_w, 0.0);
        }
        close(
            step.acoustic.incident_j - step.acoustic.outgoing_j,
            step.acoustic.deposited_j,
            2e-14,
        );
        close(
            step.acoustic.viscous_j + step.acoustic.thermal_j,
            step.acoustic.deposited_j,
            2e-14,
        );
        close(
            step.heat_report.external_heat_j,
            step.acoustic.deposited_j,
            2e-14,
        );
        close(step.heat_report.residual_j, 0.0, 2e-14);
        totals.push(step.acoustic.deposited_j);
    }
    for total in totals {
        close(
            total,
            1.3 * excitation.boundary_power().absorbed_power_w,
            2e-14,
        );
    }
}

#[test]
fn g1_native_conduction_reservoir_and_authored_heat_remain_one_thermal_balance() {
    let material = material();
    let layer = layer(&material);
    let excitation = layer
        .prepare_excitation([C64::from_re(100.0), C64::ZERO])
        .unwrap();
    let ends = [LENGTH / 2.0, LENGTH];
    let mut schedule = events(&ends, &[20.0, 0.0], &[0.3, -0.2], 1.0);
    if let Event::Birth(cell) = &mut schedule[0].event {
        cell.boundary.reservoir = Some(Reservoir {
            area_m2: AREA,
            conductance_w_k: 0.4,
            temperature_k: 299.0,
        });
    }
    let thermal = initial(schedule, &SENSIBLE);
    let heating =
        PorousHeating::try_new(&excitation, &thermal, &routes(&ends), 2, || Ok(())).unwrap();
    let step = heating
        .advance_to(&thermal, 0.2, config(1.0), || Ok(()))
        .unwrap();
    assert_eq!(step.heat_report.substeps, 1);
    // Independent explicit finite-volume balance: original thermal owner takes
    // one step here, with conduction internal and reservoir exchange external.
    let external = [0.3 + 0.4 * (299.0 - 300.1), -0.2];
    for i in 0..2 {
        let conductive = if i == 0 { -0.1 } else { 0.1 };
        let heat =
            0.2 * (heating.cells()[i].dissipation.absorbed_power_w + external[i] + conductive);
        let original = &thermal.cells()[&id(i)];
        close(
            (step.thermal.cells()[&id(i)].specific_h_j_kg - original.specific_h_j_kg)
                * original.mass_kg,
            heat,
            2e-14,
        );
        assert_eq!(step.thermal.cells()[&id(i)].boundary, original.boundary);
    }
    assert_eq!(step.thermal.links(), thermal.links());
    close(
        step.heat_report.external_heat_j,
        step.acoustic.deposited_j + 0.2 * external.iter().sum::<f64>(),
        2e-14,
    );
}

#[test]
fn g3_boundary_events_and_continuation_do_not_persist_or_duplicate_acoustic_sources() {
    let material = material();
    let layer = layer(&material);
    let excitation = layer
        .prepare_excitation([C64::from_re(100.0), C64::new(0.0, 50.0)])
        .unwrap();
    let mut schedule = events(&[LENGTH], &[0.0], &[0.1], 0.0);
    let Event::Birth(cell) = &schedule[0].event else {
        panic!("birth")
    };
    let changed = Boundary {
        source_w: 0.3,
        ..cell.boundary
    };
    schedule.push(TimedEvent {
        time_s: 0.35,
        event: Event::Boundary(id(0), changed),
    });
    let thermal = initial(schedule, &SENSIBLE);
    let heating =
        PorousHeating::try_new(&excitation, &thermal, &routes(&[LENGTH]), 1, || Ok(())).unwrap();
    let first = heating
        .advance_to(&thermal, 0.8, config(0.17), || Ok(()))
        .unwrap();
    let second = heating
        .advance_to(&first.thermal, 1.2, config(0.17), || Ok(()))
        .unwrap();
    let together = heating
        .advance_to(&thermal, 1.2, config(0.17), || Ok(()))
        .unwrap();
    let expected = 0.1 * 0.35 + 0.3 * (1.2 - 0.35) + 1.2 * heating.dissipation().absorbed_power_w;
    close(second.thermal.energy_j(), expected, 2e-14);
    close(together.thermal.energy_j(), expected, 2e-14);
    close(
        first.acoustic.deposited_j + second.acoustic.deposited_j,
        together.acoustic.deposited_j,
        2e-14,
    );
    assert_eq!(first.heat_report.events, 1);
    assert_eq!(second.thermal.cells()[&id(0)].boundary.source_w, 0.3);
    assert_eq!(thermal.cells()[&id(0)].boundary.source_w, 0.1);
}

#[test]
fn g1_latent_storage_uses_original_energy_and_domain_with_transactional_cancellation() {
    let material = material();
    let layer = layer(&material);
    let excitation = layer
        .prepare_excitation([C64::from_re(100.0), C64::ZERO])
        .unwrap();
    let thermal = initial(events(&[LENGTH], &[100.0], &[0.0], 0.0), &LATENT);
    let heating =
        PorousHeating::try_new(&excitation, &thermal, &routes(&[LENGTH]), 1, || Ok(())).unwrap();
    let latent_heat = thermal.mass_kg() * 50.0;
    let duration = latent_heat / heating.dissipation().absorbed_power_w;
    let step = heating
        .advance_to(&thermal, duration, config(duration / 7.0), || Ok(()))
        .unwrap();
    close(step.thermal.cells()[&id(0)].specific_h_j_kg, 150.0, 1e-11);
    assert_eq!(step.thermal.temperature_k(id(0)).unwrap(), 300.5);
    close(step.acoustic.deposited_j, latent_heat, 2e-13);
    assert!(matches!(
        heating.advance_to(&thermal, 5.0 * duration, config(duration / 7.0), || Ok(())),
        Err(PorousHeatingError::Thermal(DepositionError::MaterialDomain))
    ));
    let mut polls = 0;
    assert!(matches!(
        heating.advance_to(&thermal, duration, config(duration / 7.0), || {
            polls += 1;
            if polls > 40 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(PorousHeatingError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells()[&id(0)].specific_h_j_kg, 100.0);
    assert_eq!(thermal.cells()[&id(0)].boundary.source_w, 0.0);
}

#[test]
fn g0_zero_drive_preserves_native_heat_and_invalid_physical_bindings_refuse() {
    let material = material();
    let layer = layer(&material);
    let excitation = layer.prepare_excitation([C64::ZERO; 2]).unwrap();
    let ends = [LENGTH / 2.0, LENGTH];
    let schedule = events(&ends, &[20.0, 0.0], &[0.3, -0.2], 1.0);
    let thermal = initial(schedule.clone(), &SENSIBLE);
    let mapping = routes(&ends);
    let heating = PorousHeating::try_new(&excitation, &thermal, &mapping, 2, || Ok(())).unwrap();
    let step = heating
        .advance_to(&thermal, 0.8, config(0.13), || Ok(()))
        .unwrap();
    let native = thermal.advance_to(0.8, config(0.13), || Ok(())).unwrap();
    assert_eq!(step.thermal.cells(), native.0.cells());
    assert_eq!(step.heat_report, native.1);
    assert_eq!(step.acoustic.deposited_j, 0.0);
    assert_eq!(step.acoustic.residual_j, 0.0);
    let duplicated = [
        mapping[0],
        PorousHeatCell {
            cell: mapping[0].cell,
            ..mapping[1]
        },
    ];
    for bad in [&mapping[..1], &duplicated[..], &[][..]] {
        assert!(PorousHeating::try_new(&excitation, &thermal, bad, 2, || Ok(())).is_err());
    }
    assert!(PorousHeating::try_new(&excitation, &thermal, &mapping, 1, || Ok(())).is_err());
    assert!(
        PorousHeating::try_new(
            &excitation,
            &thermal,
            &mapping,
            MAX_POROUS_HEATING_CELLS + 1,
            || Ok(())
        )
        .is_err()
    );
    for mutation in [0, 1, 2, 3] {
        let mut changed = schedule.clone();
        let Event::Birth(cell) = &mut changed[1].event else {
            panic!("birth")
        };
        match mutation {
            0 => cell.geometry.volume_m3 *= 1.2,
            1 => {
                cell.geometry.start_m[2] += 0.01;
                cell.geometry.end_m[2] += 0.01;
            }
            2 => {
                cell.geometry.end_m = [
                    cell.geometry.start_m[0] - 0.8 * LENGTH / 2.0,
                    cell.geometry.start_m[1] + 0.6 * LENGTH / 2.0,
                    cell.geometry.start_m[2],
                ];
            }
            _ => cell.mass_kg *= 1.2,
        }
        let changed = initial(changed, &SENSIBLE);
        if mutation < 3 {
            assert!(PorousHeating::try_new(&excitation, &changed, &mapping, 2, || Ok(())).is_err());
        }
        assert!(
            heating
                .advance_to(&changed, 0.1, config(0.1), || Ok(()))
                .is_err()
        );
    }
    let changed_law = initial(schedule, &LATENT);
    assert!(
        heating
            .advance_to(&changed_law, 0.1, config(0.1), || Ok(()))
            .is_err()
    );
    assert!(matches!(
        PorousHeating::try_new(&excitation, &thermal, &mapping, 2, || Err(
            DepositionError::Cancelled
        )),
        Err(PorousHeatingError::Thermal(DepositionError::Cancelled))
    ));
}
