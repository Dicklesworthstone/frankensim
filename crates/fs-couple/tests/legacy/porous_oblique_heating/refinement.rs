//! G1: independent oblique field work and native thermal refinement.
//! Thermal constants are numerical fixtures, not measured material claims.

use fs_couple::porous_heating::feedback::{
    IsobaricPorousGas, PorousHarmonicDrive, PorousTemperatureCell,
    spectrum::{PorousSpectrumHeating, PorousSpectrumOptions},
};
use fs_duct::{
    compose::ComposeOptions,
    porous::{PorousLayer, PorousLayerStack},
};
use fs_material::{
    gas::{GasSpec, GasState},
    porous::{PorousCard, PorousFluidLaw},
};
use fs_math::c64::C64;
use fs_thermal::enthalpy::*;

const AREA: f64 = 0.003_906_25;
const LENGTHS: [f64; 2] = [0.031_25, 0.046_875];
const MASS: [f64; 2] = [0.0015, 0.005];
const CP: [f64; 2] = [400.0, 700.0];
const AUTHORED: [f64; 2] = [0.7, -0.15];
const CONTACT: f64 = 0.04;
const PRESSURE: f64 = 101_325.0;
const INITIAL_H: [f64; 2] = [4000.0, 42_000.0];
const IDS: [CellId; 2] = [CellId(31), CellId(7)];

struct Caloric(usize);
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h / CP[self.0];
        if t.is_finite() && (280.0..=900.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / CP[self.0]
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([91 + self.0 as u8; 32])
    }
}
static CALORIC: [Caloric; 2] = [Caloric(0), Caloric(1)];

fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-9,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:e}"
    );
}
fn thermal() -> DepositionSimulation<'static> {
    let mut events = Vec::new();
    let mut start = 0.0;
    for i in 0..2 {
        let end = start + LENGTHS[i];
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
                segment_id: i as u64,
                body_id: 23,
                material_id: MaterialId(i as u64 + 1),
                mass_kg: MASS[i],
                specific_h_j_kg: INITIAL_H[i],
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, start],
                    end_m: [0.0, 0.0, end],
                    volume_m3: AREA * LENGTHS[i],
                    surface_area_m2: 2.0 * AREA + 4.0 * AREA.sqrt() * LENGTHS[i],
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: AUTHORED[i],
                },
            }),
        });
        start = end;
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 1,
            a: IDS[0],
            b: IDS[1],
            area_m2: AREA,
            conductance_w_k: CONTACT,
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), &CALORIC[0]), (MaterialId(2), &CALORIC[1])],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 4096,
        },
    )
    .unwrap()
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0
}

fn wave(q: C64, x: f64) -> C64 {
    let phase = q.re * x;
    C64::new(phase.cos(), phase.sin()).scale((-q.im * x).exp())
}

// y = [h0,h1,outgoing_bin0_J,outgoing_bin1_J,local0_J,local1_J].
// Only the exterior S matrix comes from the native stack. Interior fields use
// independent normal-wave roots, pressure/normal-velocity continuity and plain
// exponential reconstruction. The modest attenuation makes this growing-wave
// reference safe; production instead anchors only decaying waves at each face.
// exp(-i omega t): q²=omega² rho*/K* - k_t², Z_n=omega rho*/q,
// v_t=k_t p/(omega rho*). Positive work includes BOTH velocity components.
fn rates(
    y: [f64; 6],
    materials: &[PorousFluidLaw; 2],
    spec: &GasSpec,
    drives: &[PorousHarmonicDrive; 2],
    tangential: [f64; 2],
    panels: usize,
    include_transverse_loss: bool,
) -> [f64; 6] {
    let states: [GasState; 2] =
        std::array::from_fn(|i| GasState::try_new(spec, 300.0 + y[i] / CP[i], PRESSURE).unwrap());
    let mut heat = [0.0; 2];
    let mut outgoing = [0.0; 2];
    for (bin, drive) in drives.iter().enumerate() {
        let omega = drive.angular_frequency_rad_s;
        let kt = tangential[bin];
        let layers: [PorousLayer<'_>; 2] = std::array::from_fn(|i| {
            PorousLayer::try_new_oblique(
                &materials[i],
                &states[i],
                omega,
                LENGTHS[i],
                AREA,
                drive.reference_impedance_pa_s_m3,
                kt,
            )
            .unwrap()
        });
        let stack =
            PorousLayerStack::try_new(&layers, ComposeOptions::default(), || false).unwrap();
        let s = stack.two_port().scattering;
        let a = drive.incident_pressure_pa;
        let b = [s[0] * a[0] + s[1] * a[1], s[2] * a[0] + s[3] * a[1]];
        outgoing[bin] =
            (b[0].norm_sq() + b[1].norm_sq()) / (2.0 * drive.reference_impedance_pa_s_m3);
        let mut pressure = a[0] + b[0];
        let mut velocity = (a[0] - b[0]).scale(1.0 / (drive.reference_impedance_pa_s_m3 * AREA));
        for cell in 0..2 {
            let original = materials[cell]
                .effective_fluid(&states[cell], omega)
                .unwrap();
            let rho = original.dynamic_density_kg_m3().conj();
            let bulk = original.dynamic_bulk_modulus_pa().conj();
            let q = ((rho / bulk).scale(omega * omega) - C64::from_re(kt * kt)).sqrt();
            assert!(q.re > 0.0 && q.im > 0.0);
            let zn = rho.scale(omega) / q;
            let forward = (pressure + zn * velocity).scale(0.5);
            let backward = (pressure - zn * velocity).scale(0.5);
            let viscous_weight = 0.5 * omega * rho.im;
            let thermal_weight = 0.5 * omega * bulk.recip().im;
            let dx = LENGTHS[cell] / panels as f64;
            let mut sum = 0.0;
            for index in 0..=panels {
                let x = if index == panels {
                    LENGTHS[cell]
                } else {
                    index as f64 * dx
                };
                let f = forward * wave(q, x);
                let r = backward * wave(q, -x);
                pressure = f + r;
                velocity = (f - r) / zn;
                let transverse = (pressure / rho).scale(kt / omega);
                let velocity_sq = velocity.norm_sq()
                    + if include_transverse_loss {
                        transverse.norm_sq()
                    } else {
                        0.0
                    };
                let power = viscous_weight * velocity_sq + thermal_weight * pressure.norm_sq();
                assert!(power.is_finite() && power >= 0.0);
                let weight = if index == 0 || index == panels {
                    1.0
                } else if index % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                sum += weight * power;
            }
            heat[cell] += AREA * dx * sum / 3.0;
        }
        let expected_pressure = a[1] + b[1];
        let expected_velocity =
            (b[1] - a[1]).scale(1.0 / (drive.reference_impedance_pa_s_m3 * AREA));
        assert!((pressure - expected_pressure).abs() < 2e-9 * expected_pressure.abs().max(1.0));
        assert!((velocity - expected_velocity).abs() < 2e-9 * expected_velocity.abs().max(1.0));
    }
    let contact = CONTACT * (states[1].temperature - states[0].temperature);
    [
        (heat[0] + AUTHORED[0] + contact) / MASS[0],
        (heat[1] + AUTHORED[1] - contact) / MASS[1],
        outgoing[0],
        outgoing[1],
        heat[0],
        heat[1],
    ]
}

#[test]
fn g1_oblique_field_work_and_native_heat_refine_with_transverse_viscous_loss() {
    let spec = GasSpec::dry_air_ussa1976();
    let exterior = GasState::try_new(&spec, 300.0, PRESSURE).unwrap();
    let materials = [
        PorousCard::melamine_uf_nguyen2024(),
        PorousCard::basotect_foam02(),
    ]
    .map(|card| PorousFluidLaw::try_new(card, 1.81e-5).unwrap());
    let cosines: [f64; 2] = [0.6, 0.8];
    let drives: [PorousHarmonicDrive; 2] = std::array::from_fn(|i| PorousHarmonicDrive {
        angular_frequency_rad_s: core::f64::consts::TAU * [450.0, 1100.0][i],
        area_m2: AREA,
        reference_impedance_pa_s_m3: exterior.characteristic_impedance / (AREA * cosines[i]),
        incident_pressure_pa: [
            [C64::new(600.0, -500.0), C64::new(-400.0, 900.0)],
            [C64::new(1200.0, 400.0), C64::new(600.0, -800.0)],
        ][i],
    });
    let tangential = std::array::from_fn(|i| {
        drives[i].angular_frequency_rad_s * (1.0 - cosines[i] * cosines[i]).sqrt()
            / exterior.sound_speed
    });
    let initial = thermal();
    let cells = std::array::from_fn::<_, 2, _>(|i| PorousTemperatureCell {
        cell: IDS[i],
        material: &materials[i],
    });
    let gas = IsobaricPorousGas::try_new(&spec, PRESSURE, [280.0, 900.0]).unwrap();
    let heating = PorousSpectrumHeating::try_new_oblique(
        &initial,
        &cells,
        gas,
        &[drives[1], drives[0]],
        &[tangential[1], tangential[0]],
        PorousSpectrumOptions::default(),
        || Ok(()),
    )
    .unwrap();
    let y0 = [INITIAL_H[0], INITIAL_H[1], 0.0, 0.0, 0.0, 0.0];
    let initial_rates = rates(y0, &materials, &spec, &drives, tangential, 512, true);
    let omitted = rates(y0, &materials, &spec, &drives, tangential, 512, false);
    assert!(
        (initial_rates[4] - omitted[4]) / initial_rates[4] > 0.01,
        "fixture must detect omission of transverse viscous loss"
    );
    let first_dt = 0.015625;
    let first = heating
        .advance_to(&initial, first_dt, config(first_dt), || Ok(()))
        .unwrap();
    for i in 0..2 {
        close(
            first.step.thermal.cells()[&IDS[i]].specific_h_j_kg,
            y0[i] + first_dt * initial_rates[i],
            2e-8,
        );
        close(
            first.step.cell_energy_j[i].1,
            first_dt * initial_rates[4 + i],
            3e-10,
        );
        close(
            first.bin_energy[i].outgoing_j,
            first_dt * initial_rates[2 + i],
            3e-12,
        );
    }
    let duration = 2.0;
    let mut expected = y0;
    let dt = duration / 512.0;
    for _ in 0..512 {
        let eval = |y| rates(y, &materials, &spec, &drives, tangential, 128, true);
        let a = eval(expected);
        let b = eval(std::array::from_fn(|i| expected[i] + 0.5 * dt * a[i]));
        let c = eval(std::array::from_fn(|i| expected[i] + 0.5 * dt * b[i]));
        let d = eval(std::array::from_fn(|i| expected[i] + dt * c[i]));
        expected = std::array::from_fn(|i| {
            expected[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0
        });
    }
    assert!((expected[0] - y0[0]) / CP[0] > 1.0);
    for state in [y0, expected] {
        let coarse = rates(state, &materials, &spec, &drives, tangential, 128, true);
        let fine = rates(state, &materials, &spec, &drives, tangential, 512, true);
        for i in 0..6 {
            close(coarse[i], fine[i], 5e-8 * fine[i].abs().max(1.0));
        }
    }
    let mut errors = Vec::new();
    for dt in [0.015625, 0.0078125, 0.00390625] {
        let result = heating
            .advance_to(&initial, duration, config(dt), || Ok(()))
            .unwrap();
        let mut error = [0.0; 3];
        for i in 0..2 {
            error[0] +=
                (result.step.thermal.cells()[&IDS[i]].specific_h_j_kg - expected[i]).abs() / CP[i];
            error[1] += (result.bin_energy[i].outgoing_j - expected[2 + i]).abs();
            error[2] += (result.step.cell_energy_j[i].1 - expected[4 + i]).abs();
            let incident = drives[i]
                .incident_pressure_pa
                .iter()
                .map(|p| p.norm_sq())
                .sum::<f64>()
                / (2.0 * drives[i].reference_impedance_pa_s_m3);
            close(result.bin_energy[i].incident_j, duration * incident, 3e-11);
            close(result.bin_energy[i].residual_j, 0.0, 3e-11);
            assert_eq!(
                result.heating.bins()[i].drive().angular_frequency_rad_s,
                drives[i].angular_frequency_rad_s
            );
        }
        close(
            result.step.heat_report.stored_energy_change_j,
            result.step.acoustic.deposited_j + duration * AUTHORED.iter().sum::<f64>(),
            3e-11,
        );
        close(result.step.acoustic.residual_j, 0.0, 3e-11);
        assert_eq!(
            result.work.gas_evaluations,
            2 * result.step.heat_report.substeps
        );
        assert_eq!(
            result.work.stack_solves,
            2 * result.step.heat_report.substeps
        );
        errors.push(error);
    }
    for pair in errors.windows(2) {
        for i in 0..3 {
            assert!(pair[1][i] < 0.6 * pair[0][i], "{errors:?}");
        }
    }
    assert!(
        errors[2][0] < 0.02 && errors[2][1] < 0.01 && errors[2][2] < 0.01,
        "{errors:?}"
    );
    assert_eq!(initial.time_s(), 0.0);
}
