//! G1/G3 finite heater-to-sink transport through both original counterflow ports.

use std::cell::Cell;

use fs_couple::superfluid::counterflow::{
    CounterflowHeatControls, CounterflowHeatError, CounterflowThermalContact,
    advance_counterflow_heat,
};
use fs_couple::superfluid::counterflow_pair::advance_counterflow_heat_pair;
use fs_duct::superfluid::{
    BoundaryVelocities, ChannelEnd, LinearizationLimits, MutualFriction, TwoFluidChannel,
};
use fs_material::superfluid::{NormalViscousWall, TwoFluidProperties, TwoFluidState};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalStorageBasis,
    TimedEvent,
};

const T0: f64 = 2.0;
const IDS: [CellId; 2] = [CellId(7), CellId(11)];
const MASS: [f64; 2] = [0.01, 0.02];
const CP: [f64; 2] = [5.0, 7.0];

struct Caloric {
    cp: f64,
    cubic: f64,
    lower: f64,
    upper: f64,
}
impl ThermalLaw for Caloric {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || h < self.lower || h > self.upper {
            return Err(DepositionError::MaterialDomain);
        }
        let theta = h / self.cp;
        Ok(T0 + theta + self.cubic * theta.powi(3))
    }
    fn maximum_temperature_slope(&self) -> f64 {
        let bound = self.lower.abs().max(self.upper.abs()) / self.cp;
        (1.0 + 3.0 * self.cubic * bound * bound) / self.cp
    }
}
fn laws(cubic: f64) -> [Caloric; 2] {
    CP.map(|cp| Caloric {
        cp,
        cubic,
        lower: -0.2,
        upper: 0.2,
    })
}
fn channel(theta: f64, expansion: f64) -> TwoFluidChannel {
    channel_with_wall(theta, expansion, None)
}
fn channel_with_wall(
    theta: f64,
    expansion: f64,
    wall: Option<NormalViscousWall>,
) -> TwoFluidChannel {
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
    let limits = LinearizationLimits {
        max_pressure_perturbation_pa: 1.0,
        max_temperature_perturbation_k: 0.02,
        max_component_speed_m_s: 0.1,
    };
    let mut value = match wall {
        Some(wall) => {
            TwoFluidChannel::try_new_with_normal_wall(material, 1.0, 0.01, 2, limits, wall)
        }
        None => TwoFluidChannel::try_new(material, 1.0, 0.01, 2, limits),
    }
    .unwrap();
    value
        .initialize(&[0.0; 2], &[theta; 2], &[0.0; 3], &[0.0; 3], || false)
        .unwrap();
    value
}
fn thermal<'a>(
    laws: [&'a dyn ThermalLaw; 2],
    h: [f64; 2],
    source_w: [f64; 2],
) -> DepositionSimulation<'a> {
    let events = (0..2)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64 + 1),
                mass_kg: MASS[i],
                specific_h_j_kg: h[i],
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 0.01, 0.0, 0.0],
                    volume_m3: 1e-5,
                    surface_area_m2: 1e-3,
                },
                boundary: Boundary {
                    ambient_k: T0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: source_w[i],
                },
            }),
        })
        .collect();
    let value = DepositionSimulation::new(
        events,
        &[(MaterialId(1), laws[0]), (MaterialId(2), laws[1])],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 2,
            max_substeps: 8,
        },
    )
    .unwrap();
    value
        .advance_to(0.0, controls().thermal, || Ok(()))
        .unwrap()
        .0
}
fn contacts(g: [f64; 2]) -> [CounterflowThermalContact; 2] {
    [ChannelEnd::Left, ChannelEnd::Right].map(|end| {
        let i = usize::from(end == ChannelEnd::Right);
        CounterflowThermalContact {
            end,
            cell: IDS[i],
            conductance_w_k: g[i],
        }
    })
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
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}"
    );
}
fn impedance(dt: f64) -> (f64, f64) {
    // Independent two-cell second-sound midpoint elimination, alpha=0.
    let r = 6.0 * dt * dt; // c2²=6, dx=0.5.
    let capacity = 4.0 * 16.0 * 0.01 * 0.5;
    let scale = dt / (2.0 * capacity * (1.0 + 2.0 * r));
    (scale * (1.0 + r), scale * r)
}
fn reference_heat_inventory(channel: &TwoFluidChannel) -> f64 {
    let material = channel.material();
    let volume = channel.area_m2() * channel.length_m() / channel.cells() as f64;
    (0..channel.cells())
        .map(|i| {
            let cell = channel.cell(i).unwrap();
            T0 * volume
                * (material.density_kg_m3() * cell.specific_entropy_perturbation_j_kg_k
                    + material.entropy_j_kg_k() * cell.density_perturbation_kg_m3)
        })
        .sum()
}
fn snapshot(channel: &TwoFluidChannel, thermal: &DepositionSimulation<'_>) -> Vec<f64> {
    let mut values = vec![
        channel.time_s(),
        thermal.time_s(),
        thermal.energy_j(),
        channel.cumulative_mutual_friction_heat_j(),
        channel.cumulative_normal_viscous_heat_j(),
    ];
    for id in IDS {
        values.extend([
            thermal.temperature_k(id).unwrap(),
            thermal.cells()[&id].specific_h_j_kg,
        ]);
    }
    values.extend(channel.last_mutual_friction_heat_j());
    values.extend(channel.last_normal_viscous_heat_j());
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
fn two_finite_contacts_match_independent_coupled_equations_and_availability() {
    let caloric = laws(0.0);
    for reverse_order in [false, true] {
        let theta = [0.01, -0.006];
        let channel_theta = 0.001;
        let mut wave = channel(channel_theta, 0.0);
        let mut heat = thermal(
            [&caloric[0], &caloric[1]],
            std::array::from_fn(|i| CP[i] * theta[i]),
            [0.0; 2],
        );
        let before_wave = wave.wave_energy_j().unwrap();
        let initial = reference_heat_inventory(&wave) + heat.energy_j();
        let g = [1.2, 0.7];
        let dt = 0.4 * wave.maximum_step_s();
        let (diagonal, cross) = impedance(dt);
        let a: [f64; 2] =
            std::array::from_fn(|i| 1.0 + g[i] * (diagonal + dt / (2.0 * MASS[i] * CP[i])));
        let b = g.map(|gi| gi * cross);
        let rhs: [f64; 2] = std::array::from_fn(|i| g[i] * (channel_theta - theta[i]));
        let determinant = a[0] * a[1] - b[0] * b[1];
        let expected = [
            (rhs[0] * a[1] - b[0] * rhs[1]) / determinant,
            (rhs[1] * a[0] - b[1] * rhs[0]) / determinant,
        ];
        let mut ports = contacts(g);
        if reverse_order {
            ports.swap(0, 1);
        }
        let report = advance_counterflow_heat_pair(
            &mut wave,
            &mut heat,
            ports,
            dt,
            MutualFriction::None,
            controls(),
            || false,
        )
        .unwrap();
        let mut availability = wave.wave_energy_j().unwrap() - before_wave;
        for i in 0..2 {
            near(report.heat_outward_w[i], expected[i], 2e-13);
            let after_theta = theta[i] + expected[i] * dt / (MASS[i] * CP[i]);
            near(heat.temperature_k(IDS[i]).unwrap(), T0 + after_theta, 1e-13);
            near(
                report.channel_midpoint_temperature_k[i],
                T0 + channel_theta - diagonal * expected[i] - cross * expected[1 - i],
                1e-13,
            );
            availability += 0.5 * MASS[i] * CP[i] / T0 * (after_theta.powi(2) - theta[i].powi(2));
            availability += report.contact_availability_loss_j[i];
            assert!(
                report.contact_power_residual_w[i].abs() <= controls().contact_power_tolerance_w
            );
            assert!(report.iterations[i] <= controls().max_iterations);
        }
        assert!(report.heat_outward_w[0] < 0.0 && report.heat_outward_w[1] > 0.0);
        near(availability, 0.0, 2e-15);
        near(
            reference_heat_inventory(&wave) + heat.energy_j(),
            initial,
            2e-14,
        );
        near(report.channel.boundary_mass_kg, 0.0, 1e-20);
        near(report.reference_heat_transfer_residual_j, 0.0, 1e-16);
        near(report.reference_entropy_transfer_residual_j_k, 0.0, 1e-16);
        near(report.wave_port_residual_j, 0.0, 1e-16);
    }
}

#[test]
fn repeated_nonlinear_heater_to_sink_transport_marches_native_sources_once() {
    let caloric = laws(100.0);
    for wall in [None, Some(NormalViscousWall::try_new(0.003).unwrap())] {
        let mut wave = channel_with_wall(0.0, -0.03, wall);
        let source = [0.0002, -0.0001];
        let h = [0.04, -0.05];
        let mut heat = thermal([&caloric[0], &caloric[1]], h, source);
        let initial = reference_heat_inventory(&wave) + heat.energy_j();
        let dt = 0.3 * wave.maximum_step_s();
        let mut exchange = [0.0; 2];
        let mut viscous_heat = 0.0;
        for _ in 0..80 {
            let report = advance_counterflow_heat_pair(
                &mut wave,
                &mut heat,
                contacts([0.4, 0.6]),
                dt,
                MutualFriction::None,
                controls(),
                || false,
            )
            .unwrap();
            for i in 0..2 {
                exchange[i] += report.heat_to_thermal_j[i];
                assert!(
                    report.contact_power_residual_w[i].abs()
                        <= controls().contact_power_tolerance_w
                );
            }
            near(
                report.thermal_march.external_heat_j,
                (source[0] + source[1]) * dt,
                1e-17,
            );
            near(report.reference_heat_transfer_residual_j, 0.0, 1e-16);
            near(report.reference_entropy_transfer_residual_j_k, 0.0, 1e-16);
            near(report.wave_port_residual_j, 0.0, 1e-16);
            assert_eq!(report.channel.mutual_friction_heat_j, 0.0);
            near(
                wave.last_normal_viscous_heat_j().iter().sum(),
                report.channel.normal_viscous_heat_j,
                1e-18,
            );
            viscous_heat += report.channel.normal_viscous_heat_j;
            // The original interior viscous debit remains separately exported.
            // Only the two boundary heat powers enter the finite thermal cells.
            near(
                report.heat_intake.external_heat_j,
                report.heat_to_thermal_j.iter().sum(),
                1e-16,
            );
            assert_eq!(wave.time_s(), heat.time_s());
        }
        assert!(exchange[0] < 0.0 && exchange[1] > 0.0);
        if wall.is_some() {
            assert!(viscous_heat > 0.0);
        } else {
            assert_eq!(viscous_heat, 0.0);
        }
        near(wave.cumulative_normal_viscous_heat_j(), viscous_heat, 1e-18);
        for i in 0..2 {
            near(
                MASS[i] * (heat.cells()[&IDS[i]].specific_h_j_kg - h[i]),
                exchange[i] + source[i] * heat.time_s(),
                2e-15,
            );
            assert_eq!(
                heat.temperature_k(IDS[i]).unwrap(),
                caloric[i]
                    .temperature_k(heat.cells()[&IDS[i]].specific_h_j_kg)
                    .unwrap()
            );
        }
        near(
            reference_heat_inventory(&wave) + heat.energy_j() - initial,
            (source[0] + source[1]) * heat.time_s(),
            3e-14,
        );
    }
}

#[test]
fn both_original_latent_plateaux_accept_heat_without_fitted_capacity() {
    struct Plateau;
    impl ThermalLaw for Plateau {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            if !(0.0..=0.2).contains(&h) {
                return Err(DepositionError::MaterialDomain);
            }
            Ok(T0 + (h - 0.06).max(0.0) / 5.0)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            0.2
        }
    }
    let mut wave = channel(0.005, 0.0);
    let mut heat = thermal([&Plateau, &Plateau], [0.01, 0.02], [0.0; 2]);
    let before = heat.energy_j();
    let dt = 0.2 * wave.maximum_step_s();
    let report = advance_counterflow_heat_pair(
        &mut wave,
        &mut heat,
        contacts([0.1, 0.3]),
        dt,
        MutualFriction::None,
        controls(),
        || false,
    )
    .unwrap();
    for i in 0..2 {
        assert!(report.heat_outward_w[i] > 0.0);
        assert_eq!(heat.temperature_k(IDS[i]).unwrap(), T0);
    }
    near(
        heat.energy_j() - before,
        report.heat_to_thermal_j.iter().sum(),
        1e-17,
    );
}

#[test]
fn zero_and_one_contact_reduce_to_the_original_native_equations() {
    let caloric = laws(0.0);
    for g in [[0.0, 0.0], [0.4, 0.0], [0.0, 0.4]] {
        let mut wave = channel(0.002, -0.02);
        let mut heat = thermal([&caloric[0], &caloric[1]], [0.04, -0.03], [0.001, 0.0002]);
        let mut expected_wave = channel(0.002, -0.02);
        let mut expected_heat = thermal([&caloric[0], &caloric[1]], [0.04, -0.03], [0.001, 0.0002]);
        let dt = 0.4 * wave.maximum_step_s();
        if g == [0.0; 2] {
            expected_heat = expected_heat
                .advance_to(dt, controls().thermal, || Ok(()))
                .unwrap()
                .0;
            expected_wave
                .step(dt, BoundaryVelocities::CLOSED, MutualFriction::None, || {
                    false
                })
                .unwrap();
        } else {
            let i = usize::from(g[1] > 0.0);
            advance_counterflow_heat(
                &mut expected_wave,
                &mut expected_heat,
                contacts(g)[i],
                dt,
                MutualFriction::None,
                controls(),
                || false,
            )
            .unwrap();
        }
        let report = advance_counterflow_heat_pair(
            &mut wave,
            &mut heat,
            contacts(g),
            dt,
            MutualFriction::None,
            controls(),
            || false,
        )
        .unwrap();
        for i in 0..2 {
            if g[i] == 0.0 {
                assert_eq!(report.heat_outward_w[i], 0.0);
            }
        }
        for (actual, expected) in snapshot(&wave, &heat)
            .into_iter()
            .zip(snapshot(&expected_wave, &expected_heat))
        {
            near(actual, expected, 2e-12);
        }
    }
}

#[test]
fn temporary_power_or_chart_projection_does_not_refuse_an_admissible_joint_root() {
    for chart_limit in [false, true] {
        let mut wave = channel(0.0, 0.0);
        let dt = 0.4 * wave.maximum_step_s();
        let (diagonal, cross) = impedance(dt);
        let expected_q = 0.001;
        let h: [f64; 2] = std::array::from_fn(|i| {
            -CP[i] * expected_q * (1.0 + diagonal + cross + dt / (2.0 * MASS[i] * CP[i]))
        });
        let mut caloric = laws(0.0);
        let mut c = controls();
        let limit = expected_q * 1.00001;
        if chart_limit {
            for i in 0..2 {
                caloric[i].upper = h[i] + dt * limit / MASS[i];
            }
        } else {
            c.max_abs_heat_power_w = limit;
        }
        // The isolated left-coordinate root exceeds this limit while Q_right
        // is still zero. The full positive cross response brings it inside.
        let isolated_left =
            expected_q * (1.0 + cross / (1.0 + diagonal + dt / (2.0 * MASS[0] * CP[0])));
        assert!(isolated_left > limit);
        let mut heat = thermal([&caloric[0], &caloric[1]], h, [0.0; 2]);
        let report = advance_counterflow_heat_pair(
            &mut wave,
            &mut heat,
            contacts([1.0; 2]),
            dt,
            MutualFriction::None,
            c,
            || false,
        )
        .unwrap();
        for q in report.heat_outward_w {
            near(q, expected_q, 2e-13);
        }
    }
}

#[test]
fn invalid_pair_or_failed_joint_solve_preserves_both_histories() {
    let caloric = laws(0.0);
    for case in 0..7 {
        let mut wave = channel(0.0, -0.03);
        let mut heat = thermal([&caloric[0], &caloric[1]], [0.04, -0.03], [0.0; 2]);
        let dt = 0.4 * wave.maximum_step_s();
        let mut c = controls();
        let mut ports = contacts([1.0; 2]);
        match case {
            0 => ports[1].cell = ports[0].cell,
            1 => ports[1].end = ChannelEnd::Left,
            2 => c.max_iterations = 1,
            3 => c.max_abs_heat_power_w = 1e-7,
            4 => c.max_temperature_departure_k = 1e-4,
            5 => ports[1].cell = CellId(999),
            6 => heat = heat.advance_to(dt, c.thermal, || Ok(())).unwrap().0,
            _ => unreachable!(),
        }
        let before = snapshot(&wave, &heat);
        assert!(
            advance_counterflow_heat_pair(
                &mut wave,
                &mut heat,
                ports,
                dt,
                MutualFriction::None,
                c,
                || false
            )
            .is_err()
        );
        assert_eq!(snapshot(&wave, &heat), before);
    }
}

#[test]
fn early_iterative_and_final_cancellation_are_atomic_and_retryable() {
    let caloric = laws(50.0);
    let wall = Some(NormalViscousWall::try_new(0.003).unwrap());
    let mut accepted_wave = channel_with_wall(0.0, -0.03, wall);
    let mut accepted_heat = thermal([&caloric[0], &caloric[1]], [0.04, -0.03], [0.0; 2]);
    let dt = 0.4 * accepted_wave.maximum_step_s();
    let polls = Cell::new(0);
    advance_counterflow_heat_pair(
        &mut accepted_wave,
        &mut accepted_heat,
        contacts([1.0; 2]),
        dt,
        MutualFriction::None,
        controls(),
        || {
            polls.set(polls.get() + 1);
            false
        },
    )
    .unwrap();
    for stop in [1, polls.get() / 2, polls.get()] {
        let mut wave = channel_with_wall(0.0, -0.03, wall);
        let mut heat = thermal([&caloric[0], &caloric[1]], [0.04, -0.03], [0.0; 2]);
        let before = snapshot(&wave, &heat);
        let count = Cell::new(0);
        assert_eq!(
            advance_counterflow_heat_pair(
                &mut wave,
                &mut heat,
                contacts([1.0; 2]),
                dt,
                MutualFriction::None,
                controls(),
                || {
                    count.set(count.get() + 1);
                    count.get() >= stop
                }
            ),
            Err(CounterflowHeatError::Cancelled)
        );
        assert_eq!(snapshot(&wave, &heat), before);
        advance_counterflow_heat_pair(
            &mut wave,
            &mut heat,
            contacts([1.0; 2]),
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
