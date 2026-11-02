//! G1: diffuse angular flux and native heat against independent slab fields.
//! Thermal constants are numerical fixtures, not measured material claims.

use fs_couple::porous_heating::feedback::{
    IsobaricPorousGas, PorousTemperatureCell,
    covariance::{
        PorousCrossSpectrumHeating,
        angular::{DiffuseAngularNode, DiffuseHemisphere, PorousAngularMode},
    },
    spectrum::PorousSpectrumOptions,
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
const INITIAL_H: [f64; 2] = [4000.0, 42_000.0];
const AUTHORED: [f64; 2] = [0.7, -0.15];
const CONTACT: f64 = 0.04;
const PRESSURE: f64 = 101_325.0;
const IDS: [CellId; 2] = [CellId(31), CellId(7)];
const FREQUENCIES: [f64; 2] = [450.0, 1100.0];
const IRRADIANCES: [[f64; 2]; 2] = [[1200.0, 400.0], [700.0, 900.0]];

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
        Some([131 + self.0 as u8; 32])
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
                body_id: 29,
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

#[derive(Clone, Copy)]
struct Mode {
    frequency: usize,
    mu: f64,
    weight: f64,
}
fn reference_modes(nodes: usize) -> Vec<Mode> {
    (0..2)
        .flat_map(|frequency| {
            // Increasing k_t means decreasing mu within each frequency.
            (0..nodes).rev().map(move |i| Mode {
                frequency,
                mu: (i as f64 + 0.5) / nodes as f64,
                weight: 1.0 / nodes as f64,
            })
        })
        .collect()
}
fn native_modes(exterior: &GasState, count: usize) -> Vec<PorousAngularMode> {
    let nodes: Vec<_> = (0..count)
        .map(|i| DiffuseAngularNode {
            cosine: (i as f64 + 0.5) / count as f64,
            weight_mu: 1.0 / count as f64,
        })
        .collect();
    let mut modes = Vec::new();
    // Deliberately reverse frequency and tangent ordering before admission.
    for frequency in (0..2).rev() {
        let hemisphere = DiffuseHemisphere::try_new(
            exterior,
            core::f64::consts::TAU * FREQUENCIES[frequency],
            AREA,
            IRRADIANCES[frequency],
        )
        .unwrap();
        modes.extend(hemisphere.modes(&nodes, count, || Ok(())).unwrap());
    }
    modes
}
fn wave(q: C64, x: f64) -> C64 {
    let phase = q.re * x;
    C64::new(phase.cos(), phase.sin()).scale((-q.im * x).exp())
}

// Independent left-to-right 2x2 pressure/normal-velocity transfer matrices.
// No native scattering, excitation, covariance factor or spatial-loss owner
// supplies the reference. Moderate attenuation admits this growing-wave
// calculation; production instead uses stable decaying waves at both faces.
// Result: [incident, outgoing, cell0 viscous, cell0 thermal, cell1 viscous,
// cell1 thermal], all in watts. Each hemisphere's components are independent.
fn mode_power(
    mode: Mode,
    materials: &[PorousFluidLaw; 2],
    states: &[GasState; 2],
    exterior: &GasState,
    panels: usize,
) -> [f64; 6] {
    assert_eq!(panels % 2, 0);
    let omega = core::f64::consts::TAU * FREQUENCIES[mode.frequency];
    let kt = omega * (1.0 - mode.mu * mode.mu).sqrt() / exterior.sound_speed;
    let z = exterior.characteristic_impedance / mode.mu;
    let coefficients: [(C64, C64, C64, C64); 2] = std::array::from_fn(|cell| {
        let fluid = materials[cell]
            .effective_fluid(&states[cell], omega)
            .unwrap();
        let rho = fluid.dynamic_density_kg_m3().conj();
        let bulk = fluid.dynamic_bulk_modulus_pa().conj();
        let q = ((rho / bulk).scale(omega * omega) - C64::from_re(kt * kt)).sqrt();
        assert!(q.re > 0.0 && q.im > 0.0);
        (rho, bulk, q, rho.scale(omega) / q)
    });
    let mut matrix = [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE];
    for (cell, &(_, _, q, zn)) in coefficients.iter().enumerate() {
        let forward = wave(q, LENGTHS[cell]);
        let backward = wave(q, -LENGTHS[cell]);
        let cosine = (forward + backward).scale(0.5);
        let i_sine = (forward - backward).scale(0.5);
        let m = [cosine, zn * i_sine, i_sine / zn, cosine];
        matrix = [
            m[0] * matrix[0] + m[1] * matrix[2],
            m[0] * matrix[1] + m[1] * matrix[3],
            m[2] * matrix[0] + m[3] * matrix[2],
            m[2] * matrix[1] + m[3] * matrix[3],
        ];
    }
    let mut power = [0.0; 6];
    for side in 0..2 {
        let mut incident = [C64::ZERO; 2];
        // E is hemispherical irradiance: dP = 2 A E mu dmu. Converting
        // that normal power to PEAK pressure cancels mu, with no extra 2pi.
        incident[side] = C64::from_re(
            (4.0 * exterior.characteristic_impedance
                * IRRADIANCES[mode.frequency][side]
                * mode.weight)
                .sqrt(),
        );
        let a = matrix[0] - matrix[2].scale(z);
        let b = matrix[1].scale(1.0 / z) - matrix[3];
        let reflected = (incident[1].scale(2.0) - (a + b) * incident[0]) / (a - b);
        let mut pressure = incident[0] + reflected;
        let mut velocity = (incident[0] - reflected).scale(1.0 / z);
        let transmitted = matrix[0] * pressure + matrix[1] * velocity - incident[1];
        power[0] += AREA * incident[side].norm_sq() / (2.0 * z);
        power[1] += AREA * (reflected.norm_sq() + transmitted.norm_sq()) / (2.0 * z);
        for (cell, &(rho, bulk, q, zn)) in coefficients.iter().enumerate() {
            let forward = (pressure + zn * velocity).scale(0.5);
            let backward = (pressure - zn * velocity).scale(0.5);
            let dx = LENGTHS[cell] / panels as f64;
            let mut integrals = [0.0; 2];
            for index in 0..=panels {
                let x = index as f64 * dx;
                let f = forward * wave(q, x);
                let r = backward * wave(q, -x);
                pressure = f + r;
                velocity = (f - r) / zn;
                let transverse = (pressure / rho).scale(kt / omega);
                let density = [
                    0.5 * omega * rho.im * (velocity.norm_sq() + transverse.norm_sq()),
                    0.5 * omega * bulk.recip().im * pressure.norm_sq(),
                ];
                let weight = if index == 0 || index == panels {
                    1.0
                } else if index % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                for mechanism in 0..2 {
                    assert!(density[mechanism] >= 0.0);
                    integrals[mechanism] += weight * density[mechanism];
                }
            }
            for mechanism in 0..2 {
                power[2 + 2 * cell + mechanism] += AREA * dx * integrals[mechanism] / 3.0;
            }
        }
        let exit_pressure = incident[1] + transmitted;
        let exit_velocity = (transmitted - incident[1]).scale(1.0 / z);
        assert!((pressure - exit_pressure).abs() < 2e-9 * exit_pressure.abs().max(1.0));
        assert!((velocity - exit_velocity).abs() < 2e-9 * exit_velocity.abs().max(1.0));
    }
    close(
        power[0] - power[1],
        power[2..].iter().sum(),
        2e-7 * power[0].max(1.0),
    );
    power
}

// y = [h0,h1,outgoing_frequency0_J,outgoing_frequency1_J,local0_J,local1_J].
fn rates(
    y: [f64; 6],
    materials: &[PorousFluidLaw; 2],
    spec: &GasSpec,
    exterior: &GasState,
    modes: &[Mode],
    panels: usize,
) -> [f64; 6] {
    let states =
        std::array::from_fn(|i| GasState::try_new(spec, 300.0 + y[i] / CP[i], PRESSURE).unwrap());
    let mut heat = [0.0; 2];
    let mut outgoing = [0.0; 2];
    for &mode in modes {
        let power = mode_power(mode, materials, &states, exterior, panels);
        outgoing[mode.frequency] += power[1];
        for cell in 0..2 {
            heat[cell] += power[2 + 2 * cell] + power[3 + 2 * cell];
        }
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
fn g1_independent_diffuse_modes_close_flux_local_work_and_refine_native_heat() {
    let spec = GasSpec::dry_air_ussa1976();
    let exterior = GasState::try_new(&spec, 300.0, PRESSURE).unwrap();
    let materials = [
        PorousCard::melamine_uf_nguyen2024(),
        PorousCard::basotect_foam02(),
    ]
    .map(|card| PorousFluidLaw::try_new(card, 1.81e-5).unwrap());
    let initial = thermal();
    let cells = std::array::from_fn::<_, 2, _>(|i| PorousTemperatureCell {
        cell: IDS[i],
        material: &materials[i],
    });
    let gas = IsobaricPorousGas::try_new(&spec, PRESSURE, [280.0, 900.0]).unwrap();
    let y0 = [INITIAL_H[0], INITIAL_H[1], 0.0, 0.0, 0.0, 0.0];
    let states =
        std::array::from_fn(|i| GasState::try_new(&spec, 300.0 + y0[i] / CP[i], PRESSURE).unwrap());
    let incident = AREA * IRRADIANCES.iter().flatten().sum::<f64>();
    let dense = rates(y0, &materials, &spec, &exterior, &reference_modes(256), 128);
    let mut angular_errors = Vec::new();
    for nodes in [2, 4, 8, 16] {
        let modes = native_modes(&exterior, nodes);
        let reference = reference_modes(nodes);
        let heating = PorousCrossSpectrumHeating::try_new_independent_modes(
            &initial,
            &cells,
            gas,
            &modes,
            PorousSpectrumOptions::default(),
            || Ok(()),
        )
        .unwrap();
        assert_eq!(heating.bins().len(), 2 * nodes);
        close(heating.power().incident_power_w, incident, 3e-13);
        let mut total = [0.0; 6];
        let mut outgoing = [0.0; 2];
        for (index, &mode) in reference.iter().enumerate() {
            let expected = mode_power(mode, &materials, &states, &exterior, 256);
            let bin = &heating.bins()[index];
            assert_eq!(
                bin.angular_frequency_rad_s,
                core::f64::consts::TAU * FREQUENCIES[mode.frequency]
            );
            close(
                bin.reference_impedance_pa_s_m3,
                exterior.characteristic_impedance / (AREA * mode.mu),
                3e-9,
            );
            close(
                heating.tangential_wavenumber_rad_m(index).unwrap(),
                bin.angular_frequency_rad_s * (1.0 - mode.mu * mode.mu).sqrt()
                    / exterior.sound_speed,
                3e-13,
            );
            let actual = heating.bin_power(index).unwrap();
            close(
                actual.incident_power_w,
                2.0 * AREA
                    * IRRADIANCES[mode.frequency].iter().sum::<f64>()
                    * mode.mu
                    * mode.weight,
                3e-13,
            );
            close(actual.outgoing_power_w, expected[1], 3e-11);
            for cell in 0..2 {
                let local = heating.bin_cells(index).unwrap()[cell].dissipation;
                close(local.viscous_power_w, expected[2 + 2 * cell], 3e-10);
                close(local.thermal_power_w, expected[3 + 2 * cell], 3e-10);
            }
            for i in 0..6 {
                total[i] += expected[i];
            }
            outgoing[mode.frequency] += actual.outgoing_power_w;
        }
        // Distinct angles at the SAME frequency add their independent positive
        // local powers, rather than coherently summing face pressures.
        for cell in 0..2 {
            close(
                heating.cells()[cell].dissipation.viscous_power_w,
                total[2 + 2 * cell],
                2e-9,
            );
            close(
                heating.cells()[cell].dissipation.thermal_power_w,
                total[3 + 2 * cell],
                2e-9,
            );
        }
        let angular_error = (outgoing[0] - dense[2]).abs()
            + (outgoing[1] - dense[3]).abs()
            + (total[2] + total[3] - dense[4]).abs()
            + (total[4] + total[5] - dense[5]).abs();
        angular_errors.push(angular_error);
    }
    // The two-node rule need not be in the asymptotic midpoint regime.
    assert!(angular_errors[1] < angular_errors[0], "{angular_errors:?}");
    for pair in angular_errors[1..].windows(2) {
        assert!(pair[1] < 0.4 * pair[0], "{angular_errors:?}");
    }
    assert!(angular_errors[3] < 0.001 * incident, "{angular_errors:?}");

    let modes = native_modes(&exterior, 8);
    let reference = reference_modes(8);
    let heating = PorousCrossSpectrumHeating::try_new_independent_modes(
        &initial,
        &cells,
        gas,
        &modes,
        PorousSpectrumOptions::default(),
        || Ok(()),
    )
    .unwrap();
    let first_dt = 0.015625;
    let first = heating
        .advance_to(&initial, first_dt, config(first_dt), || Ok(()))
        .unwrap();
    let initial_rates = rates(y0, &materials, &spec, &exterior, &reference, 256);
    for cell in 0..2 {
        close(
            first.step.thermal.cells()[&IDS[cell]].specific_h_j_kg,
            y0[cell] + first_dt * initial_rates[cell],
            3e-8,
        );
        close(
            first.step.cell_energy_j[cell].1,
            first_dt * initial_rates[4 + cell],
            3e-10,
        );
    }
    for (index, &mode) in reference.iter().enumerate() {
        let power = mode_power(mode, &materials, &states, &exterior, 256);
        close(
            first.bin_energy[index].incident_j,
            first_dt * power[0],
            3e-13,
        );
        close(
            first.bin_energy[index].outgoing_j,
            first_dt * power[1],
            3e-13,
        );
        close(
            first.bin_energy[index].deposited_j,
            first_dt * power[2..].iter().sum::<f64>(),
            3e-11,
        );
    }

    let duration = 2.0;
    let mut expected = y0;
    let reference_dt = duration / 256.0;
    for _ in 0..256 {
        let eval = |y| rates(y, &materials, &spec, &exterior, &reference, 64);
        let a = eval(expected);
        let b = eval(std::array::from_fn(|i| {
            expected[i] + 0.5 * reference_dt * a[i]
        }));
        let c = eval(std::array::from_fn(|i| {
            expected[i] + 0.5 * reference_dt * b[i]
        }));
        let d = eval(std::array::from_fn(|i| expected[i] + reference_dt * c[i]));
        expected = std::array::from_fn(|i| {
            expected[i] + reference_dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0
        });
    }
    assert!((expected[0] - y0[0]) / CP[0] > 1.0);
    for state in [y0, expected] {
        let coarse = rates(state, &materials, &spec, &exterior, &reference, 64);
        let fine = rates(state, &materials, &spec, &exterior, &reference, 256);
        for i in 0..6 {
            close(coarse[i], fine[i], 2e-7 * fine[i].abs().max(1.0));
        }
    }
    let mut errors = Vec::new();
    for dt in [0.03125, 0.015625, 0.0078125] {
        let result = heating
            .advance_to(&initial, duration, config(dt), || Ok(()))
            .unwrap();
        let mut error = [0.0; 3];
        let mut outgoing = [0.0; 2];
        for (index, mode) in reference.iter().enumerate() {
            outgoing[mode.frequency] += result.bin_energy[index].outgoing_j;
            close(result.bin_energy[index].residual_j, 0.0, 3e-11);
        }
        for i in 0..2 {
            error[0] +=
                (result.step.thermal.cells()[&IDS[i]].specific_h_j_kg - expected[i]).abs() / CP[i];
            error[1] += (outgoing[i] - expected[2 + i]).abs();
            error[2] += (result.step.cell_energy_j[i].1 - expected[4 + i]).abs();
        }
        close(result.step.acoustic.incident_j, duration * incident, 3e-11);
        close(result.step.acoustic.residual_j, 0.0, 3e-11);
        close(
            result.step.heat_report.stored_energy_change_j,
            result.step.acoustic.deposited_j + duration * AUTHORED.iter().sum::<f64>(),
            5e-11,
        );
        assert_eq!(
            result.work.gas_evaluations,
            2 * result.step.heat_report.substeps
        );
        assert_eq!(
            result.work.stack_solves,
            modes.len() * result.step.heat_report.substeps
        );
        errors.push(error);
    }
    for pair in errors.windows(2) {
        for i in 0..3 {
            assert!(pair[1][i] < 0.6 * pair[0][i], "{errors:?}");
        }
    }
    assert!(
        errors[2][0] < 0.03 && errors[2][1] < 0.01 && errors[2][2] < 0.01,
        "{errors:?}"
    );
    assert_eq!(initial.time_s(), 0.0);
}
