//! Connected-channel heat reaches native cells once, with joint rollback.

use std::cell::Cell;

use fs_couple::superfluid::connection::advance_connected_friction_heat;
use fs_couple::superfluid::{FrictionHeatControls, FrictionHeatError};
use fs_duct::superfluid::connection::{ConnectionDrive, ConnectionKind, TwoFluidConnection};
use fs_duct::superfluid::{
    BodyForce, CellSample, ComponentForceDensity, FaceSample, LinearizationLimits, MutualFriction,
    TwoFluidChannel,
};
use fs_material::superfluid::{NormalViscousWall, TwoFluidProperties, TwoFluidState};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalStorageBasis,
    TimedEvent,
};

const N: usize = 3;
const MASS: f64 = 0.01;
const DT: f64 = 0.004;
const ROUTES: [[CellId; N]; 2] = [
    [CellId(10), CellId(11), CellId(12)],
    [CellId(20), CellId(21), CellId(22)],
];

struct Caloric {
    upper: f64,
    positive_queries: Cell<usize>,
}
impl Caloric {
    fn new(upper: f64) -> Self {
        Self {
            upper,
            positive_queries: Cell::new(0),
        }
    }
}
impl ThermalLaw for Caloric {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h > 0.0 {
            self.positive_queries.set(self.positive_queries.get() + 1);
        }
        if !h.is_finite() || !(0.0..=self.upper).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(2.0 + h / 8.0)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.125
    }
}

fn controls() -> FrictionHeatControls {
    FrictionHeatControls {
        thermal: StepConfig {
            max_dt_s: 0.002,
            energy_tolerance_j: 1e-13,
        },
        cell_heat_tolerance_j: 1e-15,
        heat_energy_tolerance_j: 1e-14,
        combined_energy_tolerance_j: 1e-13,
    }
}

fn channels() -> (TwoFluidChannel, TwoFluidChannel, TwoFluidConnection) {
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
        max_pressure_perturbation_pa: 1.0,
        max_temperature_perturbation_k: 0.02,
        max_component_speed_m_s: 0.1,
    };
    let mut first = TwoFluidChannel::try_new(material, 1.0, 0.01, N, limits).unwrap();
    let mut second = TwoFluidChannel::try_new_with_normal_wall(
        material,
        0.7,
        0.015,
        N,
        limits,
        NormalViscousWall::try_new(0.003).unwrap(),
    )
    .unwrap();
    for (j, channel) in [&mut first, &mut second].into_iter().enumerate() {
        channel
            .initialize(
                &[0.01 * j as f64; N],
                &[0.0; N],
                &[0.0, 0.02, 0.01, 0.0],
                &[0.0, -0.007, -0.004, 0.0],
                || false,
            )
            .unwrap();
    }
    let connection = TwoFluidConnection::try_new(&first, &second, ConnectionKind::Open).unwrap();
    (first, second, connection)
}

fn drive() -> ConnectionDrive<'static> {
    ConnectionDrive {
        forces: [BodyForce::Uniform(ComponentForceDensity {
            normal_n_m3: 0.01,
            superfluid_n_m3: -0.005,
        }); 2],
        friction: [MutualFriction::Linear { gamma_kg_m3_s: 7.5 }; 2],
        ..ConnectionDrive::default()
    }
}

fn thermal(laws: &[Caloric; 2], source_w: f64) -> DepositionSimulation<'_> {
    let mut schedule = Vec::new();
    for (j, route) in ROUTES.iter().enumerate() {
        for (i, &id) in route.iter().enumerate() {
            schedule.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id,
                    segment_id: id.0,
                    body_id: j as u64,
                    material_id: MaterialId(j as u64),
                    mass_kg: MASS,
                    specific_h_j_kg: 0.0,
                    geometry: CellGeometry {
                        start_m: [i as f64, j as f64, 0.0],
                        end_m: [i as f64 + 0.1, j as f64, 0.0],
                        volume_m3: 1e-5,
                        surface_area_m2: 1e-3,
                    },
                    boundary: Boundary {
                        ambient_k: 2.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w,
                    },
                }),
            });
        }
    }
    let original = DepositionSimulation::new(
        schedule,
        &[
            (MaterialId(0), &laws[0] as &dyn ThermalLaw),
            (MaterialId(1), &laws[1] as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2 * N,
            max_links: 0,
            max_events: 2 * N,
            max_substeps: 16,
        },
    )
    .unwrap();
    original
        .advance_to(0.0, controls().thermal, || Ok(()))
        .unwrap()
        .0
}

#[derive(Debug, PartialEq)]
struct ChannelSnapshot {
    time: f64,
    energy: f64,
    momentum: f64,
    cumulative_heat: [f64; 2],
    cells: Vec<CellSample>,
    faces: Vec<FaceSample>,
    heat: Vec<f64>,
    viscous_heat: Vec<f64>,
}
fn channel_snapshot(channel: &TwoFluidChannel) -> ChannelSnapshot {
    ChannelSnapshot {
        time: channel.time_s(),
        energy: channel.wave_energy_j().unwrap(),
        momentum: channel.resolved_momentum_kg_m_s().unwrap(),
        cumulative_heat: [
            channel.cumulative_mutual_friction_heat_j(),
            channel.cumulative_normal_viscous_heat_j(),
        ],
        cells: (0..channel.cells())
            .map(|i| channel.cell(i).unwrap())
            .collect(),
        faces: (0..=channel.cells())
            .map(|i| channel.face(i).unwrap())
            .collect(),
        heat: channel.last_mutual_friction_heat_j().to_vec(),
        viscous_heat: channel.last_normal_viscous_heat_j().to_vec(),
    }
}
#[derive(Debug, PartialEq)]
struct Snapshot {
    channels: [ChannelSnapshot; 2],
    junction: [f64; 4],
    thermal_time: f64,
    thermal_cells: Vec<DepositedCell>,
}
fn snapshot(
    first: &TwoFluidChannel,
    second: &TwoFluidChannel,
    connection: &TwoFluidConnection,
    thermal: &DepositionSimulation<'_>,
) -> Snapshot {
    let flow = connection.flow_m3_s();
    Snapshot {
        channels: [channel_snapshot(first), channel_snapshot(second)],
        junction: [
            connection.time_s(),
            flow[0],
            flow[1],
            connection.kinetic_energy_j().unwrap(),
        ],
        thermal_time: thermal.time_s(),
        thermal_cells: thermal.cells().values().cloned().collect(),
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn both_connected_channels_heat_native_cells_and_sources_advance_once() {
    let laws = [Caloric::new(2.0), Caloric::new(2.0)];
    let source_w = 0.0001;
    let (mut first, mut second, mut connection) = channels();
    let mut thermal = thermal(&laws, source_w);
    for step in 1..=4 {
        let before: Vec<_> = ROUTES
            .iter()
            .flatten()
            .map(|id| thermal.cells()[id].specific_h_j_kg)
            .collect();
        let report = advance_connected_friction_heat(
            &mut first,
            &mut second,
            &mut connection,
            &mut thermal,
            [&ROUTES[0], &ROUTES[1]],
            DT,
            drive(),
            controls(),
            || false,
        )
        .unwrap();
        for (j, channel) in [&first, &second].into_iter().enumerate() {
            let wave = report.connection.channels[j];
            assert!(wave.mutual_friction_heat_j > 0.0);
            assert_eq!(report.heat[j].channel, wave);
            near(
                report.heat[j].heat_intake.external_heat_j,
                wave.dissipated_heat_j(),
                1e-14,
            );
            for (i, heat) in channel.last_dissipated_heat_j().enumerate() {
                let increase =
                    MASS * (thermal.cells()[&ROUTES[j][i]].specific_h_j_kg - before[j * N + i]);
                near(increase, source_w * DT + heat, 2e-15);
            }
        }
        assert!(report.connection.channels[1].normal_viscous_heat_j > 0.0);
        near(
            report.heat[0].thermal_march.external_heat_j,
            2.0 * N as f64 * source_w * DT,
            1e-15,
        );
        assert_eq!(report.heat[1].thermal_march.external_heat_j, 0.0);
        assert!(report.combined_storage_residual_j.abs() < 1e-13);
        assert_eq!(first.time_s(), second.time_s());
        assert_eq!(first.time_s(), connection.time_s());
        assert_eq!(first.time_s(), thermal.time_s());
        near(first.time_s(), step as f64 * DT, 1e-17);
    }
}

#[test]
fn second_caloric_refusal_preserves_all_four_accepted_owners() {
    let laws = [Caloric::new(2.0), Caloric::new(0.0)];
    let (mut first, mut second, mut connection) = channels();
    let mut thermal = thermal(&laws, 0.0);
    let before = snapshot(&first, &second, &connection, &thermal);
    let refused = advance_connected_friction_heat(
        &mut first,
        &mut second,
        &mut connection,
        &mut thermal,
        [&ROUTES[0], &ROUTES[1]],
        DT,
        drive(),
        controls(),
        || false,
    );
    assert!(matches!(
        refused,
        Err(FrictionHeatError::Thermal(DepositionError::MaterialDomain))
    ));
    assert!(laws[0].positive_queries.get() > 0);
    assert!(laws[1].positive_queries.get() > 0);
    assert_eq!(snapshot(&first, &second, &connection, &thermal), before);
}

#[test]
fn final_cancellation_preserves_both_channels_junction_and_thermal_then_retries() {
    let laws = [Caloric::new(2.0), Caloric::new(2.0)];
    let (mut reference_first, mut reference_second, mut reference_connection) = channels();
    let mut reference_thermal = thermal(&laws, 0.0001);
    let calls = Cell::new(0usize);
    advance_connected_friction_heat(
        &mut reference_first,
        &mut reference_second,
        &mut reference_connection,
        &mut reference_thermal,
        [&ROUTES[0], &ROUTES[1]],
        DT,
        drive(),
        controls(),
        || {
            calls.set(calls.get() + 1);
            false
        },
    )
    .unwrap();
    let complete_calls = calls.get();
    assert!(complete_calls > 10);
    let expected = snapshot(
        &reference_first,
        &reference_second,
        &reference_connection,
        &reference_thermal,
    );

    let (mut first, mut second, mut connection) = channels();
    let mut thermal = thermal(&laws, 0.0001);
    let before = snapshot(&first, &second, &connection, &thermal);
    calls.set(0);
    let refused = advance_connected_friction_heat(
        &mut first,
        &mut second,
        &mut connection,
        &mut thermal,
        [&ROUTES[0], &ROUTES[1]],
        DT,
        drive(),
        controls(),
        || {
            calls.set(calls.get() + 1);
            calls.get() >= complete_calls
        },
    );
    assert!(matches!(refused, Err(FrictionHeatError::Cancelled)));
    assert_eq!(calls.get(), complete_calls);
    assert_eq!(snapshot(&first, &second, &connection, &thermal), before);
    advance_connected_friction_heat(
        &mut first,
        &mut second,
        &mut connection,
        &mut thermal,
        [&ROUTES[0], &ROUTES[1]],
        DT,
        drive(),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(snapshot(&first, &second, &connection, &thermal), expected);
}

#[test]
fn shared_destinations_and_open_outer_ports_refuse_without_mutation() {
    let laws = [Caloric::new(2.0), Caloric::new(2.0)];
    let (mut first, mut second, mut connection) = channels();
    let mut thermal = thermal(&laws, 0.0);
    let before = snapshot(&first, &second, &connection, &thermal);
    let overlap = [ROUTES[0][0], ROUTES[1][1], ROUTES[1][2]];
    assert!(
        advance_connected_friction_heat(
            &mut first,
            &mut second,
            &mut connection,
            &mut thermal,
            [&ROUTES[0], &overlap],
            DT,
            drive(),
            controls(),
            || false,
        )
        .is_err()
    );
    assert_eq!(snapshot(&first, &second, &connection, &thermal), before);
    let mut open = drive();
    open.outer.left.normal_m_s = 0.001;
    assert!(
        advance_connected_friction_heat(
            &mut first,
            &mut second,
            &mut connection,
            &mut thermal,
            [&ROUTES[0], &ROUTES[1]],
            DT,
            open,
            controls(),
            || false,
        )
        .is_err()
    );
    assert_eq!(snapshot(&first, &second, &connection, &thermal), before);
}
