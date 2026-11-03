//! G1/G3/G4: original two-fluid loss enters original finite caloric cells once.

use std::cell::Cell;

use fs_couple::superfluid::{
    FrictionHeatControls, FrictionHeatError, advance_closed_friction_heat,
};
use fs_duct::superfluid::{
    CellSample, FaceSample, LinearizationLimits, MutualFriction, TwoFluidChannel,
};
use fs_material::superfluid::{TwoFluidProperties, TwoFluidState};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalStorageBasis,
    TimedEvent,
};

const N: usize = 4;
const MASS: f64 = 0.01;

struct CaloricLaw {
    heat_capacity_j_kg_k: f64,
    max_specific_energy_j_kg: f64,
}
impl ThermalLaw for CaloricLaw {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn temperature_k(&self, u: f64) -> Result<f64, DepositionError> {
        if !u.is_finite() || !(0.0..=self.max_specific_energy_j_kg).contains(&u) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(2.0 + u / self.heat_capacity_j_kg_k)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.heat_capacity_j_kg_k
    }
}

fn channel() -> TwoFluidChannel {
    let material = TwoFluidState::try_new(TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: 0.0,
        normal_viscosity_pa_s: 0.03,
    })
    .unwrap();
    let mut channel = TwoFluidChannel::try_new(
        material,
        1.0,
        0.01,
        N,
        LinearizationLimits {
            max_pressure_perturbation_pa: 1.0,
            max_temperature_perturbation_k: 0.02,
            max_component_speed_m_s: 0.1,
        },
    )
    .unwrap();
    let w: Vec<_> = (0..=N)
        .map(|i| 0.04 * (core::f64::consts::PI * i as f64 / N as f64).sin())
        .collect();
    let mut normal: Vec<_> = w.iter().map(|w| 0.75 * w).collect();
    let mut superfluid: Vec<_> = w.iter().map(|w| -0.25 * w).collect();
    normal[N] = 0.0;
    superfluid[N] = 0.0;
    channel
        .initialize(&[0.0; N], &[0.0; N], &normal, &superfluid, || false)
        .unwrap();
    channel
}

fn controls() -> FrictionHeatControls {
    FrictionHeatControls {
        thermal: StepConfig {
            max_dt_s: 0.001,
            energy_tolerance_j: 1e-13,
        },
        cell_heat_tolerance_j: 1e-15,
        heat_energy_tolerance_j: 1e-14,
        combined_energy_tolerance_j: 1e-13,
    }
}
fn destinations() -> Vec<CellId> {
    (0..N).map(|i| CellId(i as u64 + 10)).collect()
}
fn friction() -> MutualFriction {
    MutualFriction::Linear { gamma_kg_m3_s: 7.5 }
}

fn thermal(law: &dyn ThermalLaw, source_w: f64) -> DepositionSimulation<'_> {
    let schedule = destinations()
        .iter()
        .enumerate()
        .map(|(i, &id)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: i as u64,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: MASS,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [i as f64 / N as f64, 0.0, 0.0],
                    end_m: [(i + 1) as f64 / N as f64, 0.0, 0.0],
                    volume_m3: 0.01 / N as f64,
                    surface_area_m2: 0.01,
                },
                boundary: Boundary {
                    ambient_k: 2.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w,
                },
            }),
        })
        .collect();
    let original = DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: N,
            max_links: 0,
            max_events: N,
            max_substeps: 100,
        },
    )
    .unwrap();
    original
        .advance_to(0.0, controls().thermal, || Ok(()))
        .unwrap()
        .0
}

fn near(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 3e-12 * scale, "{a:e} != {b:e}");
}

#[derive(Debug, PartialEq)]
struct Snapshot {
    channel_time: f64,
    wave_energy: f64,
    cumulative_heat: f64,
    cells: Vec<CellSample>,
    faces: Vec<FaceSample>,
    heat: Vec<f64>,
    thermal_time: f64,
    thermal_cells: Vec<DepositedCell>,
}
fn snapshot(channel: &TwoFluidChannel, thermal: &DepositionSimulation<'_>) -> Snapshot {
    Snapshot {
        channel_time: channel.time_s(),
        wave_energy: channel.wave_energy_j().unwrap(),
        cumulative_heat: channel.cumulative_mutual_friction_heat_j(),
        cells: (0..channel.cells())
            .map(|i| channel.cell(i).unwrap())
            .collect(),
        faces: (0..=channel.cells())
            .map(|i| channel.face(i).unwrap())
            .collect(),
        heat: channel.last_mutual_friction_heat_j().to_vec(),
        thermal_time: thermal.time_s(),
        thermal_cells: thermal.cells().values().cloned().collect(),
    }
}

#[test]
fn physical_loss_heats_original_cells_once_across_repeated_steps() {
    // Distinct supplied calorics prove that the adapter does not copy cp from
    // the frozen acoustic reference or invent a temperature/enthalpy relation.
    let law = CaloricLaw {
        heat_capacity_j_kg_k: 8.0,
        max_specific_energy_j_kg: 2.0,
    };
    let source_w = 0.0001;
    let mut flow = channel();
    let original_material = *flow.material();
    let mut heat = thermal(&law, source_w);
    let ids = destinations();
    let dt = 0.2 * flow.maximum_step_s();
    let initial_wave_energy = flow.wave_energy_j().unwrap();
    let mut total_friction = 0.0;
    for step in 1..=2 {
        let old: Vec<_> = ids
            .iter()
            .map(|id| heat.cells()[id].specific_h_j_kg)
            .collect();
        let report = advance_closed_friction_heat(
            &mut flow,
            &mut heat,
            &ids,
            dt,
            friction(),
            controls(),
            || false,
        )
        .unwrap();
        assert!(report.channel.mutual_friction_heat_j > 0.0);
        total_friction += report.channel.mutual_friction_heat_j;
        near(
            report.thermal_march.external_heat_j,
            N as f64 * source_w * dt,
            1e-4,
        );
        near(
            report.heat_intake.external_heat_j,
            report.channel.mutual_friction_heat_j,
            1e-4,
        );
        near(report.combined_storage_residual_j, 0.0, 1e-4);
        assert_eq!(heat.time_s(), flow.time_s());
        assert!(heat.uses_material_law(MaterialId(1), &law));
        assert_eq!(
            flow.material(),
            &original_material,
            "thermal intake does not relinearize the channel"
        );
        for (i, id) in ids.iter().enumerate() {
            let u = heat.cells()[id].specific_h_j_kg;
            near(
                MASS * (u - old[i]),
                source_w * dt + flow.last_mutual_friction_heat_j()[i],
                1e-4,
            );
            near(heat.temperature_k(*id).unwrap(), 2.0 + u / 8.0, 2.0);
            assert_eq!(heat.cells()[id].boundary.source_w, source_w);
        }
        near(
            heat.energy_j(),
            step as f64 * N as f64 * source_w * dt + total_friction,
            1e-4,
        );
        near(
            initial_wave_energy - flow.wave_energy_j().unwrap(),
            total_friction,
            1e-4,
        );
        near(
            flow.cumulative_mutual_friction_heat_j(),
            total_friction,
            1e-4,
        );
    }
}

#[test]
fn thermal_domain_refusal_does_not_publish_wave_damping_or_heat() {
    let law = CaloricLaw {
        heat_capacity_j_kg_k: 16.0,
        max_specific_energy_j_kg: 0.0,
    };
    let mut flow = channel();
    let mut heat = thermal(&law, 0.0);
    let before = snapshot(&flow, &heat);
    let dt = 0.2 * flow.maximum_step_s();
    let result = advance_closed_friction_heat(
        &mut flow,
        &mut heat,
        &destinations(),
        dt,
        friction(),
        controls(),
        || false,
    );
    assert!(matches!(
        result,
        Err(FrictionHeatError::Thermal(DepositionError::MaterialDomain))
    ));
    assert_eq!(snapshot(&flow, &heat), before);
}

#[test]
fn duplicate_unknown_and_incomplete_routes_preserve_both_histories() {
    let law = CaloricLaw {
        heat_capacity_j_kg_k: 16.0,
        max_specific_energy_j_kg: 2.0,
    };
    let mut flow = channel();
    let mut heat = thermal(&law, 0.0);
    let before = snapshot(&flow, &heat);
    let dt = 0.2 * flow.maximum_step_s();
    let mut duplicate = destinations();
    duplicate[1] = duplicate[0];
    let mut unknown = destinations();
    unknown[1] = CellId(999);
    let short = vec![CellId(10)];
    for map in [duplicate, unknown, short] {
        assert!(matches!(
            advance_closed_friction_heat(
                &mut flow,
                &mut heat,
                &map,
                dt,
                friction(),
                controls(),
                || false,
            ),
            Err(FrictionHeatError::Invalid(_))
        ));
        assert_eq!(snapshot(&flow, &heat), before);
    }
    // A source owner advanced independently cannot be silently resynchronized.
    heat = heat
        .advance_to(dt, controls().thermal, || Ok(()))
        .unwrap()
        .0;
    let mismatched = snapshot(&flow, &heat);
    assert!(matches!(
        advance_closed_friction_heat(
            &mut flow,
            &mut heat,
            &destinations(),
            dt,
            friction(),
            controls(),
            || false,
        ),
        Err(FrictionHeatError::Invalid(_))
    ));
    assert_eq!(snapshot(&flow, &heat), mismatched);
}

#[test]
fn cancellation_at_entry_and_after_thermal_acceptance_rolls_back_both_owners() {
    let law = CaloricLaw {
        heat_capacity_j_kg_k: 16.0,
        max_specific_energy_j_kg: 2.0,
    };
    let mut accepted_flow = channel();
    let mut accepted_heat = thermal(&law, 0.0);
    let dt = 0.2 * accepted_flow.maximum_step_s();
    let count = Cell::new(0);
    advance_closed_friction_heat(
        &mut accepted_flow,
        &mut accepted_heat,
        &destinations(),
        dt,
        friction(),
        controls(),
        || {
            count.set(count.get() + 1);
            false
        },
    )
    .unwrap();
    let final_checkpoint = count.get();
    assert!(final_checkpoint > 5);
    for fail_at in [1, final_checkpoint] {
        let mut flow = channel();
        let mut heat = thermal(&law, 0.0);
        let before = snapshot(&flow, &heat);
        let calls = Cell::new(0);
        let result = advance_closed_friction_heat(
            &mut flow,
            &mut heat,
            &destinations(),
            dt,
            friction(),
            controls(),
            || {
                calls.set(calls.get() + 1);
                calls.get() == fail_at
            },
        );
        assert!(matches!(result, Err(FrictionHeatError::Cancelled)));
        assert_eq!(snapshot(&flow, &heat), before);
        // A retry starts from the original physical interval, not its rejected heat.
        advance_closed_friction_heat(
            &mut flow,
            &mut heat,
            &destinations(),
            dt,
            friction(),
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(
            snapshot(&flow, &heat),
            snapshot(&accepted_flow, &accepted_heat)
        );
    }
}
