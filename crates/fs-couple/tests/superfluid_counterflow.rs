//! Native finite caloric feedback drives the original two-fluid heat ports.

use std::cell::Cell;

use fs_couple::superfluid::counterflow::{
    CounterflowHeatControls, CounterflowHeatError, CounterflowThermalContact,
    advance_counterflow_heat,
};
use fs_duct::superfluid::{ChannelEnd, LinearizationLimits, MutualFriction, TwoFluidChannel};
use fs_material::superfluid::{TwoFluidProperties, TwoFluidState};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalStorageBasis,
    TimedEvent,
};

const MASS: f64 = 0.01;
const CP: f64 = 5.0;
const T0: f64 = 2.0;
const ID: CellId = CellId(7);

struct Caloric {
    cubic: f64,
    bound: f64,
}
impl ThermalLaw for Caloric {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || h.abs() > self.bound {
            return Err(DepositionError::MaterialDomain);
        }
        let theta = h / CP;
        Ok(T0 + theta + self.cubic * theta.powi(3))
    }
    fn maximum_temperature_slope(&self) -> f64 {
        (1.0 + 3.0 * self.cubic * (self.bound / CP).powi(2)) / CP
    }
}
fn channel(theta: f64, expansion: f64) -> TwoFluidChannel {
    let material = TwoFluidState::try_new(TwoFluidProperties {
        temperature_k: T0,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: expansion,
        normal_viscosity_pa_s: 0.03,
    })
    .unwrap();
    let mut value = TwoFluidChannel::try_new(
        material,
        1.0,
        0.01,
        2,
        LinearizationLimits {
            max_pressure_perturbation_pa: 1.0,
            max_temperature_perturbation_k: 0.02,
            max_component_speed_m_s: 0.1,
        },
    )
    .unwrap();
    value
        .initialize(&[0.0; 2], &[theta; 2], &[0.0; 3], &[0.0; 3], || false)
        .unwrap();
    value
}
fn thermal(law: &dyn ThermalLaw, h: f64, source_w: f64) -> DepositionSimulation<'_> {
    let value = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: ID,
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(1),
                mass_kg: MASS,
                specific_h_j_kg: h,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [0.01, 0.0, 0.0],
                    volume_m3: 1e-5,
                    surface_area_m2: 1e-3,
                },
                boundary: Boundary {
                    ambient_k: T0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w,
                },
            }),
        }],
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 8,
        },
    )
    .unwrap();
    value
        .advance_to(0.0, controls().thermal, || Ok(()))
        .unwrap()
        .0
}
fn controls() -> CounterflowHeatControls {
    CounterflowHeatControls {
        thermal: StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-13,
        },
        max_iterations: 64,
        max_abs_heat_power_w: 1.0,
        max_temperature_departure_k: 0.02,
        contact_power_tolerance_w: 1e-13,
        cell_energy_tolerance_j: 1e-15,
        transfer_energy_tolerance_j: 1e-13,
    }
}
fn contact(end: ChannelEnd, conductance: f64) -> CounterflowThermalContact {
    CounterflowThermalContact {
        end,
        cell: ID,
        conductance_w_k: conductance,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}"
    );
}
fn entropy_inventory(channel: &TwoFluidChannel) -> f64 {
    let material = channel.material();
    let volume = channel.area_m2() * channel.length_m() / channel.cells() as f64;
    (0..channel.cells())
        .map(|i| {
            let c = channel.cell(i).unwrap();
            T0 * volume
                * (material.density_kg_m3() * c.specific_entropy_perturbation_j_kg_k
                    + material.entropy_j_kg_k() * c.density_perturbation_kg_m3)
        })
        .sum()
}
fn snapshot(channel: &TwoFluidChannel, thermal: &DepositionSimulation<'_>) -> Vec<f64> {
    let mut values = vec![
        channel.time_s(),
        thermal.time_s(),
        thermal.energy_j(),
        thermal.temperature_k(ID).unwrap(),
        channel.cumulative_mutual_friction_heat_j(),
    ];
    values.extend(channel.last_mutual_friction_heat_j());
    for i in 0..channel.cells() {
        let c = channel.cell(i).unwrap();
        values.extend([
            c.pressure_perturbation_pa,
            c.temperature_perturbation_k,
            c.density_perturbation_kg_m3,
            c.specific_entropy_perturbation_j_kg_k,
        ]);
    }
    for i in 0..=channel.cells() {
        let f = channel.face(i).unwrap();
        values.extend([f.normal_velocity_m_s, f.superfluid_velocity_m_s]);
    }
    values
}

#[test]
fn finite_heat_capacity_feedback_matches_independent_two_cell_equations() {
    let law = Caloric {
        cubic: 0.0,
        bound: 0.2,
    };
    for end in [ChannelEnd::Left, ChannelEnd::Right] {
        for solid_theta in [-0.01, 0.01] {
            let mut wave = channel(0.0, 0.0);
            let mut heat = thermal(&law, CP * solid_theta, 0.0);
            let dt = 0.4 * wave.maximum_step_s();
            let conductance = 2.0;
            // Independent alpha=0 caloric midpoint elimination, including the
            // interior second-sound current and the finite solid heat capacity.
            let c2_squared = 6.0;
            let r = c2_squared * dt * dt; // dx=0.5, so dt/(2*dx)=dt.
            let fluid_cell_capacity = 4.0 * 16.0 * 0.01 * 0.5;
            let impedance = dt * (1.0 + r) / (2.0 * fluid_cell_capacity * (1.0 + 2.0 * r));
            let q = -conductance * solid_theta
                / (1.0 + conductance * (impedance + dt / (2.0 * MASS * CP)));
            let report = advance_counterflow_heat(
                &mut wave,
                &mut heat,
                contact(end, conductance),
                dt,
                MutualFriction::None,
                controls(),
                || false,
            )
            .unwrap();
            near(report.heat_outward_w, q, 2e-13);
            near(
                heat.temperature_k(ID).unwrap(),
                T0 + solid_theta + q * dt / (MASS * CP),
                1e-13,
            );
            near(
                report.channel_midpoint_temperature_k,
                T0 - impedance * q,
                1e-13,
            );
            near(report.heat_intake.stored_energy_change_j, q * dt, 1e-14);
            near(report.channel.boundary_mass_kg, 0.0, 1e-20);
            near(entropy_inventory(&wave), -q * dt, 1e-14);
            assert!(report.contact_availability_loss_j >= 0.0);
            let solid_after_theta = heat.temperature_k(ID).unwrap() - T0;
            let availability_change = wave.wave_energy_j().unwrap()
                + 0.5 * MASS * CP / T0 * (solid_after_theta.powi(2) - solid_theta.powi(2));
            near(
                availability_change + report.contact_availability_loss_j,
                0.0,
                1e-15,
            );
            assert!(
                report.channel.boundary_wave_work_j.abs() < report.heat_to_thermal_j.abs() / 100.0
            );
        }
    }
}

#[test]
fn repeated_heating_keeps_native_sources_and_reference_heat_without_replay() {
    let law = Caloric {
        cubic: 100.0,
        bound: 0.2,
    };
    let mut wave = channel(0.0, -0.05);
    let source_w = 0.0003;
    let mut heat = thermal(&law, 0.04, source_w);
    let initial = entropy_inventory(&wave) + heat.energy_j();
    let dt = 0.3 * wave.maximum_step_s();
    let mut exchange = 0.0;
    for _ in 0..60 {
        let report = advance_counterflow_heat(
            &mut wave,
            &mut heat,
            contact(ChannelEnd::Left, 0.3),
            dt,
            MutualFriction::None,
            controls(),
            || false,
        )
        .unwrap();
        exchange += report.heat_to_thermal_j;
        near(report.thermal_march.external_heat_j, source_w * dt, 1e-17);
        near(report.reference_heat_transfer_residual_j, 0.0, 1e-17);
        near(report.reference_entropy_transfer_residual_j_k, 0.0, 1e-17);
        near(report.wave_port_residual_j, 0.0, 1e-17);
        assert!(report.contact_power_residual_w.abs() <= controls().contact_power_tolerance_w);
        assert_eq!(wave.time_s(), heat.time_s());
    }
    near(
        heat.energy_j() - MASS * 0.04,
        exchange + source_w * heat.time_s(),
        2e-15,
    );
    near(
        entropy_inventory(&wave) + heat.energy_j() - initial,
        source_w * heat.time_s(),
        2e-14,
    );
    let h = heat.cells()[&ID].specific_h_j_kg;
    assert_eq!(
        heat.temperature_k(ID).unwrap(),
        law.temperature_k(h).unwrap()
    );
    assert!(wave.cell(0).unwrap().temperature_perturbation_k.abs() > 1e-5);
}

#[test]
fn original_latent_plateau_accepts_heat_without_an_invented_heat_capacity() {
    struct Plateau;
    impl ThermalLaw for Plateau {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            if !(0.0..=0.2).contains(&h) {
                return Err(DepositionError::MaterialDomain);
            }
            Ok(T0 + (h - 0.05).max(0.0) / CP)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            1.0 / CP
        }
    }
    let mut wave = channel(0.005, 0.0);
    let mut heat = thermal(&Plateau, 0.01, 0.0);
    let dt = 0.2 * wave.maximum_step_s();
    let before = heat.energy_j();
    let report = advance_counterflow_heat(
        &mut wave,
        &mut heat,
        contact(ChannelEnd::Right, 0.1),
        dt,
        MutualFriction::None,
        controls(),
        || false,
    )
    .unwrap();
    assert!(report.heat_outward_w > 0.0);
    assert!(heat.energy_j() > before);
    assert_eq!(heat.temperature_k(ID).unwrap(), T0);
    near(heat.energy_j() - before, report.heat_to_thermal_j, 1e-17);
}

#[test]
fn zero_conductance_is_closed_wave_and_preserves_ordinary_thermal_sources() {
    let law = Caloric {
        cubic: 0.0,
        bound: 0.2,
    };
    let mut wave = channel(0.002, -0.02);
    let mut heat = thermal(&law, 0.04, 0.001);
    let mut expected_wave = channel(0.002, -0.02);
    let dt = 0.4 * wave.maximum_step_s();
    let expected_heat = heat
        .advance_to(dt, controls().thermal, || Ok(()))
        .unwrap()
        .0;
    expected_wave
        .step(
            dt,
            fs_duct::superfluid::BoundaryVelocities::CLOSED,
            MutualFriction::None,
            || false,
        )
        .unwrap();
    let report = advance_counterflow_heat(
        &mut wave,
        &mut heat,
        contact(ChannelEnd::Left, 0.0),
        dt,
        MutualFriction::None,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(report.heat_to_thermal_j, 0.0);
    assert_eq!(
        snapshot(&wave, &heat),
        snapshot(&expected_wave, &expected_heat)
    );
}

#[test]
fn power_solve_domain_and_clock_refusals_leave_both_histories_unchanged() {
    let law = Caloric {
        cubic: 0.0,
        bound: 0.2,
    };
    for case in 0..5 {
        let mut wave = channel(0.0, -0.03);
        let mut heat = thermal(&law, 0.04, 0.0);
        let dt = 0.4 * wave.maximum_step_s();
        let mut c = controls();
        let mut port = contact(ChannelEnd::Left, 1.0);
        match case {
            0 => c.max_abs_heat_power_w = 1e-6,
            1 => c.max_iterations = 1,
            2 => c.max_temperature_departure_k = 1e-4,
            3 => port.cell = CellId(999),
            4 => {
                heat = heat.advance_to(dt, c.thermal, || Ok(())).unwrap().0;
            }
            _ => unreachable!(),
        }
        let before = snapshot(&wave, &heat);
        assert!(
            advance_counterflow_heat(
                &mut wave,
                &mut heat,
                port,
                dt,
                MutualFriction::None,
                c,
                || false
            )
            .is_err()
        );
        assert_eq!(snapshot(&wave, &heat), before);
    }
    // A caloric chart that admits only the original source state cannot give
    // up heat. The bracket must stop at that chart rather than extrapolate it.
    struct Point;
    impl ThermalLaw for Point {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            if h == 0.0 {
                Ok(2.01)
            } else {
                Err(DepositionError::MaterialDomain)
            }
        }
        fn maximum_temperature_slope(&self) -> f64 {
            0.0
        }
    }
    let mut wave = channel(0.0, 0.0);
    let mut heat = thermal(&Point, 0.0, 0.0);
    let dt = 0.4 * wave.maximum_step_s();
    let before = snapshot(&wave, &heat);
    assert!(matches!(
        advance_counterflow_heat(
            &mut wave,
            &mut heat,
            contact(ChannelEnd::Left, 1.0),
            dt,
            MutualFriction::None,
            controls(),
            || false
        ),
        Err(CounterflowHeatError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
    assert_eq!(snapshot(&wave, &heat), before);
}

#[test]
fn early_mid_solve_and_final_cancellation_are_atomic_and_retryable() {
    let law = Caloric {
        cubic: 10.0,
        bound: 0.2,
    };
    let port = contact(ChannelEnd::Left, 1.0);
    let mut accepted_wave = channel(0.0, -0.03);
    let mut accepted_heat = thermal(&law, 0.04, 0.0);
    let dt = 0.4 * accepted_wave.maximum_step_s();
    let count = Cell::new(0);
    advance_counterflow_heat(
        &mut accepted_wave,
        &mut accepted_heat,
        port,
        dt,
        MutualFriction::None,
        controls(),
        || {
            count.set(count.get() + 1);
            false
        },
    )
    .unwrap();
    for stop in [1, count.get() / 2, count.get()] {
        let mut wave = channel(0.0, -0.03);
        let mut heat = thermal(&law, 0.04, 0.0);
        let before = snapshot(&wave, &heat);
        let polls = Cell::new(0);
        assert_eq!(
            advance_counterflow_heat(
                &mut wave,
                &mut heat,
                port,
                dt,
                MutualFriction::None,
                controls(),
                || {
                    polls.set(polls.get() + 1);
                    polls.get() >= stop
                }
            ),
            Err(CounterflowHeatError::Cancelled)
        );
        assert_eq!(snapshot(&wave, &heat), before);
        advance_counterflow_heat(
            &mut wave,
            &mut heat,
            port,
            dt,
            MutualFriction::None,
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(
            snapshot(&wave, &heat),
            snapshot(&accepted_wave, &accepted_heat)
        );
    }
}
