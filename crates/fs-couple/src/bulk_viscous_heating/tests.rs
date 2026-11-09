//! Authored constitutive fixtures establish coupling, not material calibration.
#[path = "finite_window/tests.rs"]
mod finite_window;

use super::*;
use fs_duct::{compose::ComposeOptions, fluid::bulk::BulkTractionCylinder};
use fs_material::{
    fluid_acoustics::{
        FluidAcousticState,
        thermoviscous::{ThermoviscousFluidState, absorption::BulkThermoviscousEquilibrium},
    },
    fluid_transport::{FluidTransportState, ThermalFluidTransportState},
};
use fs_math::c64::C64;
use fs_thermal::enthalpy::{
    Boundary, DepositedCell, DepositionBudget, ThermalLaw, ThermalLink, TimedEvent,
};

const T0: f64 = 12.5;
const RHO: f64 = 2.0;
const RADIUS: f64 = 0.02;
fn area() -> f64 {
    core::f64::consts::PI * RADIUS * RADIUS
}
fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-12,
    }
}
fn near(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * b.abs().max(1e-280),
        "{a:.17e} != {b:.17e}"
    );
}

struct Caloric {
    cv: f64,
    kind: u8,
    identity: Option<[u8; 32]>,
    basis: ThermalStorageBasis,
}
impl Caloric {
    fn new(cv: f64, kind: u8) -> Self {
        Self {
            cv,
            kind,
            identity: Some([kind; 32]),
            basis: ThermalStorageBasis::InternalEnergy,
        }
    }
}
impl ThermalLaw for Caloric {
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.basis
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        self.identity
    }
    fn temperature_k(&self, u: f64) -> Result<f64, DepositionError> {
        if !u.is_finite() || u < 0.0 || u > 10.0 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(T0
            + match self.kind {
                // Positive nonlinear Cv = cv + (T-T0), integrated before inversion.
                1 => 2.0 * u / ((self.cv * self.cv + 2.0 * u).sqrt() + self.cv),
                // Original constant-temperature latent interval 1 <= u <= 3.
                2 => {
                    if u < 1.0 {
                        u / self.cv
                    } else if u <= 3.0 {
                        1.0 / self.cv
                    } else {
                        (u - 2.0) / self.cv
                    }
                }
                _ => u / self.cv,
            })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.cv
    }
}
fn response(
    alpha: f64,
    length: f64,
    omega: f64,
    drive: C64,
    matched: bool,
) -> BulkTractionAbsorptionResponse {
    let transport = FluidTransportState::try_new(RHO, 8.0, 0.6, 0.003).unwrap();
    let beta_s = 1.0 / 200.0;
    let beta_t = beta_s + alpha * alpha * T0 / 16.0;
    let medium = BulkThermoviscousEquilibrium::from_states(
        ThermoviscousFluidState::try_new(
            FluidAcousticState::try_new(RHO, 10.0).unwrap(),
            transport,
            8.0 * beta_s / beta_t,
        )
        .unwrap(),
        ThermalFluidTransportState::try_new(transport, T0, Some(alpha)).unwrap(),
        0.002,
    )
    .unwrap();
    let cylinder = BulkTractionCylinder::new(medium.absorption(), RADIUS, length, omega).unwrap();
    let load = matched.then(|| {
        cylinder
            .wave()
            .traction_impedance_pa_s_m()
            .scale(1.0 / area())
    });
    cylinder
        .absorption_response(&medium, load, drive, ComposeOptions::default())
        .unwrap()
}
fn id(i: usize) -> CellId {
    CellId(80 - 3 * i as u64)
}
fn point(x: f64) -> [f64; 3] {
    [0.6 * x, 0.8 * x, 0.0]
}
fn routes(ends: &[f64]) -> Vec<BulkHeatCell> {
    ends.iter()
        .enumerate()
        .map(|(i, &end_distance_m)| BulkHeatCell {
            cell: id(i),
            end_distance_m,
        })
        .collect()
}
fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: T0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    }
}
fn events(ends: &[f64], conduction: f64) -> Vec<TimedEvent> {
    let mut result = Vec::new();
    let mut start = 0.0;
    for (i, &end) in ends.iter().enumerate() {
        let volume = area() * (end - start);
        result.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: id(i),
                segment_id: i as u64,
                body_id: 7,
                material_id: MaterialId(9),
                mass_kg: RHO * volume,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: point(start),
                    end_m: point(end),
                    volume_m3: volume,
                    surface_area_m2: 2.0 * area()
                        + 2.0 * core::f64::consts::PI * RADIUS * (end - start),
                },
                boundary: boundary(0.0),
            }),
        });
        if i > 0 && conduction > 0.0 {
            result.push(TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: i as u64,
                    a: id(i - 1),
                    b: id(i),
                    area_m2: area(),
                    conductance_w_k: conduction,
                }),
            });
        }
        start = end;
    }
    result
}
fn initial(events: Vec<TimedEvent>, law: &dyn ThermalLaw) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        events,
        &[(MaterialId(9), law)],
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
fn prepare<'a>(
    response: &'a BulkTractionAbsorptionResponse,
    thermal: &DepositionSimulation<'_>,
    ends: &[f64],
) -> BulkViscousHeating<'a> {
    BulkViscousHeating::try_new(response, thermal, &routes(ends), 16, 1.0, || Ok(())).unwrap()
}

// Independent stress reconstruction: t = p - eta du/dx. Integrate its actual
// continuum viscous heat, without the production interval/Gram operator.
fn viscous_quadrature(response: &BulkTractionAbsorptionResponse, start: f64, end: f64) -> f64 {
    let panels = 2048;
    let dx = (end - start) / panels as f64;
    let eta = 4.0 * 0.003 / 3.0 + 0.002;
    let mut sum = 0.0;
    for i in 0..=panels {
        let field = response.at_distance(start + i as f64 * dx).unwrap();
        let heat = (field.pressure_pa - field.normal_traction_pa).norm_sq() / (2.0 * eta);
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight * heat;
    }
    sum * dx * area() / 3.0
}

#[test]
fn coherent_local_heat_reaches_original_cells_once_on_the_native_clock() {
    let law = Caloric::new(8.0, 0);
    let wave = response(0.0, 0.04, 1000.0, C64::new(2e-5, -1e-5), false);
    let ends = [0.009, 0.027, 0.04];
    let start = initial(events(&ends, 0.0), &law);
    let source = prepare(&wave, &start, &ends);
    assert!(core::ptr::eq(source.response(), &wave));
    for cell in source.cells() {
        near(
            cell.absorption.viscous_heat_w,
            viscous_quadrature(&wave, cell.start_distance_m, cell.end_distance_m),
            4e-12,
        );
    }
    let duration = 0.031;
    let step = source
        .advance_to(&start, duration, config(0.004), || Ok(()))
        .unwrap();
    let mut total = 0.0;
    for ((id, energy), cell) in step.cell_energy_j.iter().zip(source.cells()) {
        near(*energy, cell.absorption.viscous_heat_w * duration, 2e-15);
        near(
            step.thermal.cells()[id].mass_kg * step.thermal.cells()[id].specific_h_j_kg,
            *energy,
            3e-15,
        );
        assert_eq!(step.thermal.cells()[id].boundary.source_w, 0.0);
        total += energy;
    }
    near(step.heat_report.external_heat_j, total, 3e-15);
    near(step.heat_report.stored_energy_change_j, total, 3e-15);
    near(step.acoustic.viscous_heat_j, total, 2e-15);
    assert_eq!(step.acoustic.conductive_entropy_j_k, 0.0);
    assert!(step.acoustic.balance_residual_j.abs() < 2e-12 * total);
    let midway = source
        .advance_to(&start, 0.012, config(0.004), || Ok(()))
        .unwrap();
    let resumed = source
        .advance_to(&midway.thermal, duration, config(0.004), || Ok(()))
        .unwrap();
    near(resumed.thermal.energy_j(), step.thermal.energy_j(), 3e-15);
    let one_start = initial(events(&[0.04], 0.0), &law);
    let one = prepare(&wave, &one_start, &[0.04])
        .advance_to(&one_start, duration, config(0.004), || Ok(()))
        .unwrap();
    near(one.acoustic.viscous_heat_j, total, 3e-15);
}

#[test]
fn conductive_availability_is_not_deposited_as_caloric_heat() {
    let law = Caloric::new(4.0, 0);
    let wave = response(0.08, 0.01, 10_000.0, C64::new(2e-5, -1e-5), true);
    let start = initial(events(&[0.003, 0.01], 0.0), &law);
    let source = prepare(&wave, &start, &[0.003, 0.01]);
    let step = source
        .advance_to(&start, 0.02, config(0.003), || Ok(()))
        .unwrap();
    assert!(step.acoustic.thermal_availability_loss_j > 100.0 * step.acoustic.viscous_heat_j);
    near(
        step.acoustic.thermal_availability_loss_j,
        T0 * step.acoustic.conductive_entropy_j_k,
        2e-15,
    );
    near(
        step.heat_report.external_heat_j,
        step.acoustic.viscous_heat_j,
        2e-15,
    );
    near(step.thermal.energy_j(), step.acoustic.viscous_heat_j, 2e-15);
    let flux_drop = step.acoustic.left_acoustic_work_j - step.acoustic.right_acoustic_work_j;
    near(
        flux_drop,
        step.acoustic.viscous_heat_j + step.acoustic.thermal_availability_loss_j,
        3e-13,
    );
    for (id, joules) in step.cell_energy_j {
        let cell = &step.thermal.cells()[&id];
        near(
            step.thermal.temperature_k(id).unwrap(),
            T0 + joules / (cell.mass_kg * law.cv),
            2e-15,
        );
    }
}

#[test]
fn native_conduction_and_scheduled_sources_share_one_thermal_account() {
    let law = Caloric::new(8.0, 0);
    let wave = response(0.0, 0.04, 1000.0, C64::new(2e-5, -1e-5), false);
    let ends = [0.02, 0.04];
    let mut schedule = events(&ends, 0.01);
    schedule.push(TimedEvent {
        time_s: 0.013,
        event: Event::Boundary(id(0), boundary(2e-5)),
    });
    let start = initial(schedule, &law);
    let step = prepare(&wave, &start, &ends)
        .advance_to(&start, 0.04, config(0.008), || Ok(()))
        .unwrap();
    let authored = 2e-5 * (0.04 - 0.013);
    near(
        step.heat_report.external_heat_j,
        step.acoustic.viscous_heat_j + authored,
        3e-15,
    );
    near(
        step.thermal.energy_j(),
        step.heat_report.external_heat_j,
        3e-15,
    );
    assert_eq!(step.heat_report.events, 1);
    assert_eq!(step.thermal.cells()[&id(0)].boundary.source_w, 2e-5);
    assert_eq!(step.thermal.links()[&1].conductance_w_k, 0.01);
    let mut isolated_schedule = events(&ends, 0.0);
    isolated_schedule.push(TimedEvent {
        time_s: 0.013,
        event: Event::Boundary(id(0), boundary(2e-5)),
    });
    let isolated_start = initial(isolated_schedule, &law);
    let isolated = prepare(&wave, &isolated_start, &ends)
        .advance_to(&isolated_start, 0.04, config(0.008), || Ok(()))
        .unwrap();
    near(isolated.thermal.energy_j(), step.thermal.energy_j(), 3e-15);
    let contrast = |s: &DepositionSimulation<'_>| {
        (s.temperature_k(id(0)).unwrap() - s.temperature_k(id(1)).unwrap()).abs()
    };
    assert!(contrast(&step.thermal) < contrast(&isolated.thermal));
}

#[test]
fn original_nonlinear_and_latent_internal_energy_laws_remain_authoritative() {
    let wave = response(0.0, 0.04, 1000.0, C64::new(2e-5, -1e-5), false);
    for kind in [1, 2] {
        let law = Caloric::new(8.0, kind);
        let start = initial(events(&[0.04], 0.0), &law);
        let source = prepare(&wave, &start, &[0.04]);
        let mass = start.cells()[&id(0)].mass_kg;
        let duration = mass * 2.0 / source.cells()[0].absorption.viscous_heat_w;
        let step = source
            .advance_to(&start, duration, config(duration / 8.0), || Ok(()))
            .unwrap();
        near(step.thermal.cells()[&id(0)].specific_h_j_kg, 2.0, 2e-15);
        near(
            step.thermal.temperature_k(id(0)).unwrap(),
            law.temperature_k(2.0).unwrap(),
            2e-15,
        );
        near(step.heat_report.stored_energy_change_j, mass * 2.0, 2e-15);
        if kind == 2 {
            assert_eq!(step.thermal.temperature_k(id(0)).unwrap(), T0 + 1.0 / 8.0);
        }
    }
}

#[test]
fn mismatched_volume_mass_temperature_caloric_basis_and_identity_refuse() {
    let wave = response(0.0, 0.04, 1000.0, C64::new(2e-5, -1e-5), false);
    let law = Caloric::new(8.0, 0);
    let admit = |thermal: &DepositionSimulation<'_>| {
        BulkViscousHeating::try_new(&wave, thermal, &routes(&[0.04]), 4, 1.0, || Ok(())).is_ok()
    };
    for change in 0..4 {
        let mut schedule = events(&[0.04], 0.0);
        let Event::Birth(cell) = &mut schedule[0].event else {
            unreachable!()
        };
        match change {
            0 => cell.mass_kg *= 1.01,
            1 => cell.geometry.volume_m3 *= 1.01,
            2 => cell.geometry.end_m[0] += 0.001,
            _ => cell.specific_h_j_kg = 0.1,
        }
        assert!(!admit(&initial(schedule, &law)));
    }
    let mut wrong = Caloric::new(8.0, 0);
    wrong.basis = ThermalStorageBasis::Enthalpy;
    assert!(!admit(&initial(events(&[0.04], 0.0), &wrong)));
    wrong.basis = ThermalStorageBasis::InternalEnergy;
    wrong.identity = None;
    assert!(!admit(&initial(events(&[0.04], 0.0), &wrong)));
    let start = initial(events(&[0.04], 0.0), &law);
    let source = prepare(&wave, &start, &[0.04]);
    wrong.identity = Some([77; 32]);
    let changed = initial(events(&[0.04], 0.0), &wrong);
    assert!(
        source
            .advance_to(&changed, 0.02, config(0.002), || Ok(()))
            .is_err()
    );
    for partition in [
        routes(&[0.03]),
        vec![
            BulkHeatCell {
                cell: id(0),
                end_distance_m: 0.02,
            },
            BulkHeatCell {
                cell: id(0),
                end_distance_m: 0.04,
            },
        ],
    ] {
        assert!(BulkViscousHeating::try_new(&wave, &start, &partition, 4, 1.0, || Ok(())).is_err());
    }
}

#[test]
fn temperature_domain_and_cancellation_refusals_publish_no_candidate() {
    let law = Caloric::new(8.0, 0);
    let wave = response(0.0, 0.04, 1000.0, C64::new(2e-5, -1e-5), false);
    let start = initial(events(&[0.04], 0.0), &law);
    let source =
        BulkViscousHeating::try_new(&wave, &start, &routes(&[0.04]), 4, 1e-8, || Ok(())).unwrap();
    assert!(matches!(
        source.advance_to(&start, 0.02, config(0.002), || Ok(())),
        Err(BulkViscousHeatingError::TemperatureExcursion { .. })
    ));
    let source = prepare(&wave, &start, &[0.04]);
    let mut checkpoints = 0;
    assert!(matches!(
        source.advance_to(&start, 0.02, config(0.002), || {
            checkpoints += 1;
            if checkpoints == 25 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(BulkViscousHeatingError::Thermal(DepositionError::Cancelled))
    ));
    let power = source.cells()[0].absorption.viscous_heat_w;
    let outside_domain = start.cells()[&id(0)].mass_kg * 20.0 / power;
    assert!(matches!(
        source.advance_to(&start, outside_domain, config(outside_domain), || Ok(())),
        Err(BulkViscousHeatingError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
    assert_eq!(start.time_s(), 0.0);
    assert_eq!(start.energy_j(), 0.0);
    assert_eq!(start.cells()[&id(0)].boundary.source_w, 0.0);
}

#[test]
fn silent_drive_deposits_exact_zero_without_changing_original_thermal_sources() {
    let law = Caloric::new(8.0, 0);
    let wave = response(0.08, 0.01, 10_000.0, C64::ZERO, true);
    let start = initial(events(&[0.01], 0.0), &law);
    let step = prepare(&wave, &start, &[0.01])
        .advance_to(&start, 0.02, config(0.003), || Ok(()))
        .unwrap();
    assert_eq!(step.thermal.energy_j(), 0.0);
    assert_eq!(step.acoustic.viscous_heat_j, 0.0);
    assert_eq!(step.acoustic.conductive_entropy_j_k, 0.0);
    assert_eq!(step.acoustic.balance_residual_j, 0.0);
}
