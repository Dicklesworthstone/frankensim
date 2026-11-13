//! Direct covariance-field oracle and shared native thermal refinement.
//! The authored thermal constants are numerical fixtures, not material claims.

use fs_couple::porous_heating::feedback::{IsobaricPorousGas, PorousTemperatureCell};
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
        Some([61 + self.0 as u8; 32])
    }
}
static CALORIC: [Caloric; 2] = [Caloric(0), Caloric(1)];

fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1.0e-9,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:.3e}"
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
                body_id: 20,
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

use fs_couple::porous_heating::feedback::covariance::{
    PorousCovarianceBin, PorousCrossSpectrumHeating,
};
use fs_couple::porous_heating::feedback::spectrum::PorousSpectrumOptions;
use fs_duct::porous::PressureCovariance;

const FREQUENCIES_HZ: [f64; 2] = [450.0, 1100.0];

fn matrices() -> [[C64; 4]; 2] {
    [
        [
            C64::from_re(1.5e6),
            C64::new(0.45e6, 0.35e6),
            C64::new(0.45e6, -0.35e6),
            C64::from_re(0.8e6),
        ],
        [
            C64::from_re(0.4e6),
            C64::new(-0.18e6, 0.13e6),
            C64::new(-0.18e6, -0.13e6),
            C64::from_re(0.6e6),
        ],
    ]
}

// The authored matrices are well inside the positive-semidefinite cone. This
// direct Hermitian quadratic form deliberately never factors or diagonalizes C.
// C_LR = E[p_L p_R*]; H is the row mapping incoming pressure to one output field.
fn quadratic(h: [C64; 2], c: [C64; 4]) -> f64 {
    let value =
        c[0].re * h[0].norm_sq() + c[3].re * h[1].norm_sq() + 2.0 * (h[0] * c[1] * h[1].conj()).re;
    assert!(value.is_finite() && value >= 0.0);
    value
}

// y = [h0, h1, outgoing_bin0_J, outgoing_bin1_J, local0_J, local1_J].
// Original unit left/right excitations provide H for pressure and superficial
// velocity. Independent Simpson integration of H C H* supplies local work.
// C is RMS, so the native peak-field factor 1/2 cancels the conversion 2*C:
// q = -omega Im(rho) E[|v|²] - omega Im(1/K) E[|p|²].
fn derivative(
    y: [f64; 6],
    materials: &[PorousFluidLaw; 2],
    spec: &GasSpec,
    covariance: &[[C64; 4]; 2],
    reference_impedance: f64,
    panels: usize,
) -> [f64; 6] {
    let states: [GasState; 2] =
        std::array::from_fn(|i| GasState::try_new(spec, 300.0 + y[i] / CP[i], PRESSURE).unwrap());
    let mut heat = [0.0; 2];
    let mut outgoing = [0.0; 2];
    for bin in 0..2 {
        let omega = core::f64::consts::TAU * FREQUENCIES_HZ[bin];
        let layers: [PorousLayer<'_>; 2] = std::array::from_fn(|i| {
            PorousLayer::try_new(
                &materials[i],
                &states[i],
                omega,
                LENGTHS[i],
                AREA,
                reference_impedance,
            )
            .unwrap()
        });
        let stack =
            PorousLayerStack::try_new(&layers, ComposeOptions::default(), || false).unwrap();
        let left = stack
            .prepare_excitation([C64::ONE, C64::ZERO], || false)
            .unwrap();
        let right = stack
            .prepare_excitation([C64::ZERO, C64::ONE], || false)
            .unwrap();
        let s = stack.two_port().scattering;
        outgoing[bin] = (quadratic([s[0], s[1]], covariance[bin])
            + quadratic([s[2], s[3]], covariance[bin]))
            / reference_impedance;
        for cell in 0..2 {
            let response = layers[cell].material_response();
            let viscous = -omega * response.dynamic_density_kg_m3().im;
            let thermal = -omega * response.dynamic_bulk_modulus_pa().recip().im;
            let dx = LENGTHS[cell] / panels as f64;
            let mut sum = 0.0;
            for index in 0..=panels {
                let x = if index == panels {
                    LENGTHS[cell]
                } else {
                    index as f64 * dx
                };
                let a = left.layer_excitations()[cell].at_depth(x).unwrap();
                let b = right.layer_excitations()[cell].at_depth(x).unwrap();
                let power = viscous
                    * quadratic(
                        [a.superficial_velocity_m_s, b.superficial_velocity_m_s],
                        covariance[bin],
                    )
                    + thermal * quadratic([a.pressure_pa, b.pressure_pa], covariance[bin]);
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

fn reference(
    materials: &[PorousFluidLaw; 2],
    spec: &GasSpec,
    covariance: &[[C64; 4]; 2],
    impedance: f64,
    duration: f64,
) -> [f64; 6] {
    let mut y = [INITIAL_H[0], INITIAL_H[1], 0.0, 0.0, 0.0, 0.0];
    let dt = duration / 512.0;
    for _ in 0..512 {
        let k1 = derivative(y, materials, spec, covariance, impedance, 128);
        let k2 = derivative(
            std::array::from_fn(|i| y[i] + 0.5 * dt * k1[i]),
            materials,
            spec,
            covariance,
            impedance,
            128,
        );
        let k3 = derivative(
            std::array::from_fn(|i| y[i] + 0.5 * dt * k2[i]),
            materials,
            spec,
            covariance,
            impedance,
            128,
        );
        let k4 = derivative(
            std::array::from_fn(|i| y[i] + dt * k3[i]),
            materials,
            spec,
            covariance,
            impedance,
            128,
        );
        y = std::array::from_fn(|i| y[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    y
}

#[test]
fn g1_partial_coherence_matches_direct_field_covariance_and_shared_native_thermal_refinement() {
    let spec = GasSpec::dry_air_ussa1976();
    let exterior = GasState::try_new(&spec, 300.0, PRESSURE).unwrap();
    let impedance = exterior.characteristic_impedance / AREA;
    let materials = [
        PorousCard::melamine_uf_nguyen2024(),
        PorousCard::basotect_foam02(),
    ]
    .map(|card| PorousFluidLaw::try_new(card, 1.81e-5).unwrap());
    let covariance = matrices();
    let mut bins: [PorousCovarianceBin; 2] = std::array::from_fn(|i| PorousCovarianceBin {
        angular_frequency_rad_s: core::f64::consts::TAU * FREQUENCIES_HZ[i],
        area_m2: AREA,
        reference_impedance_pa_s_m3: impedance,
        covariance: PressureCovariance::try_new(covariance[i]).unwrap(),
    });
    bins.reverse();
    let initial = thermal();
    let cells = std::array::from_fn::<_, 2, _>(|i| PorousTemperatureCell {
        cell: IDS[i],
        material: &materials[i],
    });
    let gas = IsobaricPorousGas::try_new(&spec, PRESSURE, [280.0, 900.0]).unwrap();
    let heating = PorousCrossSpectrumHeating::try_new(
        &initial,
        &cells,
        gas,
        &bins,
        PorousSpectrumOptions::default(),
        || Ok(()),
    )
    .unwrap();
    let y0 = [INITIAL_H[0], INITIAL_H[1], 0.0, 0.0, 0.0, 0.0];
    let rates0 = derivative(y0, &materials, &spec, &covariance, impedance, 512);
    let reversed = covariance.map(|c| [c[0], c[2], c[1], c[3]]);
    let wrong_phase = derivative(y0, &materials, &spec, &reversed, impedance, 128);
    assert!(
        (wrong_phase[4] - rates0[4]).abs() > 1e-3,
        "fixture must distinguish C_LR from its conjugate"
    );
    let dt = 0.015625;
    let first = heating
        .advance_to(&initial, dt, config(dt), || Ok(()))
        .unwrap();
    let y1 = std::array::from_fn(|i| y0[i] + dt * rates0[i]);
    let rates1 = derivative(y1, &materials, &spec, &covariance, impedance, 512);
    let second = first
        .heating
        .advance_to(&first.step.thermal, 2.0 * dt, config(dt), || Ok(()))
        .unwrap();
    for i in 0..2 {
        close(
            first.step.thermal.cells()[&IDS[i]].specific_h_j_kg,
            y1[i],
            1e-8,
        );
        close(
            second.step.thermal.cells()[&IDS[i]].specific_h_j_kg,
            y1[i] + dt * rates1[i],
            2e-8,
        );
        close(first.step.cell_energy_j[i].1, dt * rates0[4 + i], 2e-11);
        close(second.step.cell_energy_j[i].1, dt * rates1[4 + i], 2e-11);
        close(first.bin_energy[i].outgoing_j, dt * rates0[2 + i], 2e-13);
        close(second.bin_energy[i].outgoing_j, dt * rates1[2 + i], 2e-11);
    }
    let duration = 4.0;
    let expected = reference(&materials, &spec, &covariance, impedance, duration);
    assert!((expected[0] - INITIAL_H[0]) / CP[0] > 20.0);
    // Bound the independent spatial-oracle error well below the tested native
    // time errors at both ends of this smooth, moderately attenuating fixture.
    for y in [y0, expected] {
        let coarse = derivative(y, &materials, &spec, &covariance, impedance, 128);
        let fine = derivative(y, &materials, &spec, &covariance, impedance, 512);
        for i in 0..6 {
            close(coarse[i], fine[i], 3e-8 * fine[i].abs().max(1.0));
        }
    }
    let mut errors = Vec::new();
    for dt in [0.015625, 0.0078125, 0.00390625] {
        let result = heating
            .advance_to(&initial, duration, config(dt), || Ok(()))
            .unwrap();
        let step = &result.step;
        let mut error = [0.0; 3];
        for i in 0..2 {
            error[0] += (step.thermal.cells()[&IDS[i]].specific_h_j_kg - expected[i]).abs() / CP[i];
            error[1] += (result.bin_energy[i].outgoing_j - expected[2 + i]).abs();
            error[2] += (step.cell_energy_j[i].1 - expected[4 + i]).abs();
            let incoming = (covariance[i][0].re + covariance[i][3].re) / impedance;
            close(result.bin_energy[i].incident_j, duration * incoming, 3e-11);
            close(result.bin_energy[i].residual_j, 0.0, 3e-11);
            assert_eq!(
                result.heating.bins()[i].angular_frequency_rad_s,
                core::f64::consts::TAU * FREQUENCIES_HZ[i]
            );
        }
        close(
            step.acoustic.outgoing_j,
            result
                .bin_energy
                .iter()
                .map(|energy| energy.outgoing_j)
                .sum::<f64>(),
            3e-11,
        );
        close(
            step.acoustic.deposited_j,
            result
                .bin_energy
                .iter()
                .map(|energy| energy.deposited_j)
                .sum::<f64>(),
            3e-11,
        );
        let authored = duration * AUTHORED.iter().sum::<f64>();
        close(
            step.heat_report.external_heat_j,
            step.acoustic.deposited_j + authored,
            3e-11,
        );
        close(
            step.heat_report.stored_energy_change_j,
            step.acoustic.deposited_j + authored,
            3e-11,
        );
        close(step.acoustic.residual_j, 0.0, 3e-11);
        assert_eq!(step.heat_report.substeps, (duration / dt) as usize);
        assert_eq!(result.work.gas_evaluations, 2 * step.heat_report.substeps);
        assert_eq!(result.work.stack_solves, 2 * step.heat_report.substeps);
        errors.push(error);
    }
    for pair in errors.windows(2) {
        for component in 0..3 {
            assert!(
                pair[1][component] < 0.6 * pair[0][component],
                "partial-coherence feedback must converge: {errors:?}"
            );
        }
    }
    assert!(
        errors[2][0] < 0.02,
        "aggregate temperature error [K]: {errors:?}"
    );
    assert!(
        errors[2][1] < 0.01 && errors[2][2] < 0.01,
        "outgoing and local work error [J]: {errors:?}"
    );
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(heating.gas_states()[0].temperature, 310.0);
}
