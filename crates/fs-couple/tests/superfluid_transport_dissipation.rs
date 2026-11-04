//! Heat-driven superfluid transport with explicitly mapped distributed losses.

use fs_couple::superfluid::counterflow::{
    CounterflowHeatControls, CounterflowHeatError, CounterflowThermalContact,
};
use fs_couple::superfluid::counterflow_pair::dissipation::{
    CounterflowDissipationStep, advance_counterflow_heat_pair_with_dissipation,
};
use fs_duct::superfluid::{ChannelEnd, LinearizationLimits, MutualFriction, TwoFluidChannel};
use fs_material::superfluid::{
    NormalViscousSpectrum, NormalViscousWall, TwoFluidProperties, TwoFluidState,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalStorageBasis,
    TimedEvent,
};
use std::cell::Cell;

const IDS: [CellId; 4] = [CellId(7), CellId(11), CellId(101), CellId(102)];
const MASS: [f64; 4] = [0.01, 0.02, 0.01, 0.015];
const CP: [f64; 4] = [5.0, 7.0, 4.0, 6.0];
const DT: f64 = 0.005;

struct Caloric {
    cp: f64,
    upper: f64,
}
impl ThermalLaw for Caloric {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || h < -1.0 || h > self.upper {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(2.0 + h / self.cp)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.cp
    }
}

fn channel(profile: bool) -> TwoFluidChannel {
    let material = TwoFluidState::try_new(TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: -0.1,
        normal_viscosity_pa_s: 0.03,
    })
    .unwrap();
    let limits = LinearizationLimits {
        max_pressure_perturbation_pa: 5.0,
        max_temperature_perturbation_k: 0.02,
        max_component_speed_m_s: 0.1,
    };
    let mut value = if profile {
        TwoFluidChannel::try_new_with_normal_profile(
            material,
            1.0,
            0.01,
            2,
            limits,
            NormalViscousSpectrum::parallel_plates(material, 0.2, 4).unwrap(),
        )
    } else {
        TwoFluidChannel::try_new_with_normal_wall(
            material,
            1.0,
            0.01,
            2,
            limits,
            NormalViscousWall::try_new(0.003).unwrap(),
        )
    }
    .unwrap();
    value
        .initialize(
            &[0.0; 2],
            &[0.0; 2],
            &[0.0, 0.02, 0.0],
            &[0.0, -0.006, 0.0],
            || false,
        )
        .unwrap();
    value
}
fn controls() -> CounterflowHeatControls {
    CounterflowHeatControls {
        thermal: StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-13,
        },
        max_iterations: 96,
        max_abs_heat_power_w: 0.01,
        max_temperature_departure_k: 0.02,
        contact_power_tolerance_w: 1e-13,
        cell_energy_tolerance_j: 1e-15,
        transfer_energy_tolerance_j: 1e-13,
    }
}
fn thermal<'a>(laws: &'a [Caloric; 4], sources: [f64; 4]) -> DepositionSimulation<'a> {
    let h = [0.01, -0.014, 0.0, 0.0];
    let events = (0..4)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: MASS[i],
                specific_h_j_kg: h[i],
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 0.1, 0.0, 0.0],
                    volume_m3: 1e-5,
                    surface_area_m2: 1e-3,
                },
                boundary: Boundary {
                    ambient_k: 2.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: sources[i],
                },
            }),
        })
        .collect();
    let references: Vec<_> = laws
        .iter()
        .enumerate()
        .map(|(i, law)| (MaterialId(i as u64), law as &dyn ThermalLaw))
        .collect();
    DepositionSimulation::new(
        events,
        &references,
        DepositionBudget {
            max_cells: 4,
            max_links: 0,
            max_events: 4,
            max_substeps: 4,
        },
    )
    .unwrap()
    .advance_to(0.0, controls().thermal, || Ok(()))
    .unwrap()
    .0
}
fn contacts() -> [CounterflowThermalContact; 2] {
    [ChannelEnd::Left, ChannelEnd::Right].map(|end| CounterflowThermalContact {
        end,
        cell: IDS[usize::from(end == ChannelEnd::Right)],
        conductance_w_k: 0.2,
    })
}
fn advance(
    wave: &mut TwoFluidChannel,
    heat: &mut DepositionSimulation<'_>,
    route: &[CellId],
    cancelled: impl Fn() -> bool,
) -> Result<CounterflowDissipationStep, CounterflowHeatError> {
    advance_counterflow_heat_pair_with_dissipation(
        wave,
        heat,
        contacts(),
        route,
        DT,
        MutualFriction::Linear { gamma_kg_m3_s: 3.0 },
        controls(),
        cancelled,
    )
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tol={tolerance:e}"
    );
}
fn snapshot(channel: &TwoFluidChannel, thermal: &DepositionSimulation<'_>) -> Vec<f64> {
    let mut result = vec![
        channel.time_s(),
        channel.wave_energy_j().unwrap(),
        channel.normal_profile_excess_energy_j().unwrap(),
        channel.cumulative_mutual_friction_heat_j(),
        channel.cumulative_normal_viscous_heat_j(),
        thermal.time_s(),
        thermal.energy_j(),
    ];
    for id in IDS {
        result.extend([
            thermal.cells()[&id].specific_h_j_kg,
            thermal.temperature_k(id).unwrap(),
        ]);
    }
    for i in 0..channel.cells() {
        let c = channel.cell(i).unwrap();
        result.extend([
            c.pressure_perturbation_pa,
            c.temperature_perturbation_k,
            c.density_perturbation_kg_m3,
            c.specific_entropy_perturbation_j_kg_k,
        ]);
    }
    for i in 0..=channel.cells() {
        let f = channel.face(i).unwrap();
        result.extend([f.normal_velocity_m_s, f.superfluid_velocity_m_s]);
        let mut m = 0;
        while let Some(value) = channel.normal_profile_mode_velocity_m_s(i, m) {
            result.push(value);
            m += 1;
        }
    }
    result.extend(channel.last_mutual_friction_heat_j());
    result.extend(channel.last_normal_viscous_heat_j());
    result
}

#[test]
fn finite_contacts_and_distributed_native_calorics_receive_each_transfer_once() {
    let laws = CP.map(|cp| Caloric { cp, upper: 1.0 });
    let sources = [1e-5, -2e-5, 1e-6, 2e-6];
    for profile in [false, true] {
        let mut wave = channel(profile);
        let mut heat = thermal(&laws, sources);
        let mut total_loss = 0.0;
        for _ in 0..20 {
            let h = IDS.map(|id| heat.cells()[&id].specific_h_j_kg);
            let temperatures = IDS.map(|id| heat.temperature_k(id).unwrap());
            let before_energy = heat.energy_j();
            let report = advance(&mut wave, &mut heat, &IDS[2..], || false).unwrap();
            let pair = report.contacts;
            let spatial: Vec<_> = wave.last_dissipated_heat_j().collect();
            let loss: f64 = spatial.iter().sum();
            assert!(pair.channel.normal_viscous_heat_j > 0.0);
            assert!(pair.channel.mutual_friction_heat_j > 0.0);
            for i in 0..4 {
                let intake = if i < 2 {
                    pair.heat_to_thermal_j[i]
                } else {
                    spatial[i - 2]
                };
                let actual = MASS[i] * (heat.cells()[&IDS[i]].specific_h_j_kg - h[i]);
                near(actual, sources[i] * DT + intake, 2e-15);
                if i < 2 {
                    let midpoint = 0.5 * (temperatures[i] + heat.temperature_k(IDS[i]).unwrap());
                    near(pair.thermal_midpoint_temperature_k[i], midpoint, 2e-14);
                    near(
                        pair.heat_outward_w[i],
                        0.2 * (pair.channel_midpoint_temperature_k[i] - midpoint),
                        2e-13,
                    );
                }
            }
            let contact_heat: f64 = pair.heat_to_thermal_j.iter().sum();
            near(
                pair.thermal_march.external_heat_j,
                sources.iter().sum::<f64>() * DT,
                2e-15,
            );
            near(pair.heat_intake.external_heat_j, contact_heat, 2e-15);
            near(report.dissipation_intake.external_heat_j, loss, 2e-15);
            near(report.dissipated_heat_j, loss, 2e-15);
            near(
                heat.energy_j() - before_energy,
                sources.iter().sum::<f64>() * DT + contact_heat + loss,
                2e-15,
            );
            assert!(report.dissipation_transfer_residual_j.abs() <= 1e-13);
            assert!(pair.reference_heat_transfer_residual_j.abs() <= 1e-13);
            total_loss += loss;
        }
        near(
            total_loss,
            wave.cumulative_mutual_friction_heat_j() + wave.cumulative_normal_viscous_heat_j(),
            2e-15,
        );
        if profile {
            assert!(wave.normal_profile_excess_energy_j().unwrap() > 0.0);
        }
    }
}

#[test]
fn aliases_incomplete_routes_and_receiving_caloric_refusal_preserve_all_owners() {
    let laws = CP.map(|cp| Caloric { cp, upper: 1.0 });
    let mut wave = channel(true);
    let mut heat = thermal(&laws, [0.0; 4]);
    advance(&mut wave, &mut heat, &IDS[2..], || false).unwrap();
    let accepted = snapshot(&wave, &heat);
    for route in [
        vec![IDS[0], IDS[3]],
        vec![IDS[2], IDS[2]],
        vec![IDS[2]],
        vec![IDS[2], CellId(999)],
    ] {
        assert!(advance(&mut wave, &mut heat, &route, || false).is_err());
        assert_eq!(snapshot(&wave, &heat), accepted);
    }
    let mut refusing = CP.map(|cp| Caloric { cp, upper: 1.0 });
    refusing[2].upper = 0.0;
    let mut wave = channel(true);
    let mut heat = thermal(&refusing, [0.0; 4]);
    let accepted = snapshot(&wave, &heat);
    assert!(matches!(
        advance(&mut wave, &mut heat, &IDS[2..], || false),
        Err(CounterflowHeatError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
    assert_eq!(snapshot(&wave, &heat), accepted);
}

#[test]
fn early_and_final_cancellation_preserve_profile_memory_and_both_heat_intakes() {
    let laws = CP.map(|cp| Caloric { cp, upper: 1.0 });
    let mut reference_wave = channel(true);
    let mut reference_heat = thermal(&laws, [0.0; 4]);
    advance(&mut reference_wave, &mut reference_heat, &IDS[2..], || {
        false
    })
    .unwrap();
    assert!(reference_wave.normal_profile_excess_energy_j().unwrap() > 0.0);
    let count = Cell::new(0_usize);
    let expected = advance(&mut reference_wave, &mut reference_heat, &IDS[2..], || {
        count.set(count.get() + 1);
        false
    })
    .unwrap();
    let reference = snapshot(&reference_wave, &reference_heat);
    assert!(count.get() > 8);
    for at in [1, count.get() / 2, count.get()] {
        let mut wave = channel(true);
        let mut heat = thermal(&laws, [0.0; 4]);
        advance(&mut wave, &mut heat, &IDS[2..], || false).unwrap();
        let accepted = snapshot(&wave, &heat);
        let polls = Cell::new(0_usize);
        let failed = advance(&mut wave, &mut heat, &IDS[2..], || {
            polls.set(polls.get() + 1);
            polls.get() >= at
        });
        assert!(failed.is_err(), "cancellation {at} was not admitted");
        assert_eq!(snapshot(&wave, &heat), accepted);
        let retried = advance(&mut wave, &mut heat, &IDS[2..], || false).unwrap();
        assert_eq!(retried, expected);
        assert_eq!(snapshot(&wave, &heat), reference);
    }
}
