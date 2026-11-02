//! Time-refinement oracle using original material and acoustic owners directly.
//! The authored thermal constants are numerical fixtures, not material claims.

use fs_couple::porous_heating::feedback::{
    IsobaricPorousGas, PorousHarmonicDrive, PorousTemperatureCell, TemperaturePorousHeating,
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

// y = [h_left, h_right, outgoing_J, absorbed_left_J, absorbed_right_J].
// This bypasses TemperaturePorousHeating, its gas/packet caches, source maps,
// heat accounting and native stepper. Only the ORIGINAL constitutive and
// coherent acoustic owners supply power; thermal balance is authored here.
fn derivative(
    y: [f64; 5],
    materials: &[PorousFluidLaw; 2],
    spec: &GasSpec,
    drive: PorousHarmonicDrive,
) -> [f64; 5] {
    let temperatures = [300.0 + y[0] / CP[0], 300.0 + y[1] / CP[1]];
    let layers: [PorousLayer<'_>; 2] = std::array::from_fn(|i| {
        let gas = GasState::try_new(spec, temperatures[i], PRESSURE).unwrap();
        PorousLayer::try_new(
            &materials[i],
            &gas,
            drive.angular_frequency_rad_s,
            LENGTHS[i],
            AREA,
            drive.reference_impedance_pa_s_m3,
        )
        .unwrap()
    });
    let stack = PorousLayerStack::try_new(&layers, ComposeOptions::default(), || false).unwrap();
    let field = stack
        .prepare_excitation(drive.incident_pressure_pa, || false)
        .unwrap();
    let heat: [f64; 2] = std::array::from_fn(|i| {
        field.layer_excitations()[i]
            .dissipation_between(0.0, LENGTHS[i])
            .unwrap()
            .absorbed_power_w
    });
    let exchange = CONTACT * (temperatures[1] - temperatures[0]);
    [
        (heat[0] + AUTHORED[0] + exchange) / MASS[0],
        (heat[1] + AUTHORED[1] - exchange) / MASS[1],
        field.boundary_power().outgoing_power_w,
        heat[0],
        heat[1],
    ]
}

fn reference(
    materials: &[PorousFluidLaw; 2],
    spec: &GasSpec,
    drive: PorousHarmonicDrive,
    duration: f64,
) -> [f64; 5] {
    let mut y = [INITIAL_H[0], INITIAL_H[1], 0.0, 0.0, 0.0];
    let dt = duration / 1024.0;
    for _ in 0..1024 {
        let k1 = derivative(y, materials, spec, drive);
        let k2 = derivative(
            std::array::from_fn(|i| y[i] + 0.5 * dt * k1[i]),
            materials,
            spec,
            drive,
        );
        let k3 = derivative(
            std::array::from_fn(|i| y[i] + 0.5 * dt * k2[i]),
            materials,
            spec,
            drive,
        );
        let k4 = derivative(
            std::array::from_fn(|i| y[i] + dt * k3[i]),
            materials,
            spec,
            drive,
        );
        y = std::array::from_fn(|i| y[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    y
}

#[test]
fn g1_native_porous_feedback_and_outgoing_energy_refine_to_independent_coherent_ode() {
    let spec = GasSpec::dry_air_ussa1976();
    let exterior = GasState::try_new(&spec, 300.0, PRESSURE).unwrap();
    let materials = [
        PorousCard::melamine_uf_nguyen2024(),
        PorousCard::basotect_foam02(),
    ]
    .map(|card| PorousFluidLaw::try_new(card, 1.81e-5).unwrap());
    let drive = PorousHarmonicDrive {
        angular_frequency_rad_s: core::f64::consts::TAU * 1100.0,
        area_m2: AREA,
        reference_impedance_pa_s_m3: exterior.characteristic_impedance / AREA,
        incident_pressure_pa: [C64::new(1200.0, 400.0), C64::new(600.0, -800.0)],
    };
    let initial = thermal();
    let gas = IsobaricPorousGas::try_new(&spec, PRESSURE, [280.0, 900.0]).unwrap();
    let cells = std::array::from_fn::<_, 2, _>(|i| PorousTemperatureCell {
        cell: IDS[i],
        material: &materials[i],
    });
    let heating = TemperaturePorousHeating::try_new(
        &initial,
        &cells,
        gas,
        drive,
        ComposeOptions::default(),
        || Ok(()),
    )
    .unwrap();
    let y0 = [INITIAL_H[0], INITIAL_H[1], 0.0, 0.0, 0.0];
    let rates0 = derivative(y0, &materials, &spec, drive);
    let dt = 0.03125;
    let first = heating
        .advance_to(&initial, dt, config(dt), || Ok(()))
        .unwrap();
    let y1 = std::array::from_fn(|i| y0[i] + dt * rates0[i]);
    let rates1 = derivative(y1, &materials, &spec, drive);
    let second = first
        .heating
        .advance_to(&first.step.thermal, 2.0 * dt, config(dt), || Ok(()))
        .unwrap();
    for i in 0..2 {
        close(
            first.step.thermal.cells()[&IDS[i]].specific_h_j_kg,
            y1[i],
            3e-10,
        );
        close(
            second.step.thermal.cells()[&IDS[i]].specific_h_j_kg,
            y1[i] + dt * rates1[i],
            3e-10,
        );
        close(first.step.cell_energy_j[i].1, dt * rates0[3 + i], 3e-13);
        close(second.step.cell_energy_j[i].1, dt * rates1[3 + i], 3e-13);
    }
    close(first.step.acoustic.outgoing_j, dt * rates0[2], 3e-13);
    close(second.step.acoustic.outgoing_j, dt * rates1[2], 3e-13);

    let duration = 8.0;
    let expected = reference(&materials, &spec, drive, duration);
    let final_rates = derivative(expected, &materials, &spec, drive);
    assert!((expected[0] - INITIAL_H[0]) / CP[0] > 20.0);
    assert!(
        (final_rates[3] - rates0[3]).abs() > 0.01,
        "fixture must exercise changing local acoustic loss"
    );
    let mut errors = Vec::new();
    for dt in [0.0625, 0.03125, 0.015625] {
        let result = heating
            .advance_to(&initial, duration, config(dt), || Ok(()))
            .unwrap();
        let step = &result.step;
        let mut temperature_error = 0.0;
        let mut acoustic_error = (step.acoustic.outgoing_j - expected[2]).abs();
        for i in 0..2 {
            let after = &step.thermal.cells()[&IDS[i]];
            temperature_error += (after.specific_h_j_kg - expected[i]).abs() / CP[i];
            acoustic_error += (step.cell_energy_j[i].1 - expected[3 + i]).abs();
            assert_eq!(after.mass_kg, MASS[i]);
            assert_eq!(after.material_id, MaterialId(i as u64 + 1));
            assert_eq!(after.boundary.source_w, AUTHORED[i]);
        }
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
        close(
            step.acoustic.incident_j,
            duration * heating.boundary_power().incident_power_w,
            3e-11,
        );
        assert_eq!(step.heat_report.substeps, (duration / dt) as usize);
        assert_eq!(result.work.stack_solves, step.heat_report.substeps);
        assert_eq!(result.work.gas_evaluations, 2 * step.heat_report.substeps);
        errors.push([temperature_error, acoustic_error]);
    }
    for pair in errors.windows(2) {
        for component in 0..2 {
            assert!(
                pair[1][component] < 0.6 * pair[0][component],
                "accepted-start feedback must converge with refinement: {errors:?}"
            );
        }
    }
    assert!(
        errors[2][0] < 0.02,
        "aggregate temperature error [K]: {errors:?}"
    );
    assert!(
        errors[2][1] < 0.01,
        "aggregate acoustic work error [J]: {errors:?}"
    );
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(heating.gas_states()[0].temperature, 310.0);
}
