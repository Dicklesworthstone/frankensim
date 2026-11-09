//! G1/G3 pressure-driven flow through the existing MAC momentum and heat owners.

use super::*;
use crate::conjugate::{
    EnergyConfig, ThermalFace, ThermalSetup, UnsteadyConfig, simple_unsteady,
    simple_unsteady_with_pressure_schedule, solve_energy,
};

fn fluid() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 2.0,
        conductivity_w_m_k: 0.01,
        kinematic_viscosity_m2_s: 0.25,
    }
}

fn channel(n: usize, axis: usize, inlet: f64, outlet: f64) -> (VoxelDomain, SimpleConfig) {
    let transverse = (axis + 1) % 3;
    let mut widths = [vec![1.0], vec![1.0], vec![1.0]];
    widths[axis] = vec![0.25; 8];
    widths[transverse] = vec![1.0 / n as f64; n];
    let domain = VoxelDomain::graded(widths).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[2 * axis] = FvBoundary::Pressure { pressure_pa: inlet };
    faces[2 * axis + 1] = FvBoundary::Pressure {
        pressure_pa: outlet,
    };
    faces[2 * transverse] = FvBoundary::wall();
    faces[2 * transverse + 1] = FvBoundary::wall();
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-10;
    config.pressure_tolerance = 1e-12;
    (domain, config)
}

#[test]
fn pressure_opening_poiseuille_refines_and_preserves_the_face_pressure_datum() {
    let gate = CancelGate::new();
    let mut errors = Vec::new();
    for n in [4, 8] {
        let (domain, config) = channel(n, 0, 7.5, 7.0);
        let flow = simple_flow(&domain, &fluid(), &config, &gate).unwrap();
        // Fully developed planar flow, H=1, L=2, mu=1/4, depth=1:
        // u(y)=dp/(2 mu L)*y*(H-y), Q=dp*H^3/(12 mu L).
        let exact_q = 0.5 / (12.0 * 0.25 * 2.0);
        let q = flow.field.boundary_outflow(&domain, Face3::XMax);
        errors.push((q - exact_q).abs());
        assert!(q > 0.0);
        assert!((q / exact_q - (1.0 + 2.0 / (n * n) as f64)).abs() < 1e-7);
        assert!((flow.report.inflow_m3_s - q).abs() < 1e-12);
        assert!((flow.report.outflow_m3_s - q).abs() < 1e-12);
        assert!(flow.report.max_divergence_m3_s < 1e-12);
        for c in 0..domain.cell_count() {
            let at = domain.coords(c);
            let x = domain.center(at[0], at[1], at[2])[0];
            assert!((flow.pressure_pa[c] - (7.5 - 0.25 * x)).abs() < 2e-7);
        }
    }
    assert!((errors[0] / errors[1] - 4.0).abs() < 1e-5, "{errors:?}");
}

#[test]
fn pressure_opening_offsets_and_equal_pressure_do_not_create_motion() {
    let gate = CancelGate::new();
    let (domain, config) = channel(4, 0, 0.5, 0.0);
    let base = simple_flow(&domain, &fluid(), &config, &gate).unwrap();
    for offset in [-1024.0, 1_048_576.0] {
        let (_, shifted) = channel(4, 0, 0.5 + offset, offset);
        let flow = simple_flow(&domain, &fluid(), &shifted, &gate).unwrap();
        assert_eq!(flow.velocity_m_s, base.velocity_m_s);
        assert_eq!(flow.report, base.report);
        for (pressure, original) in flow.pressure_pa.iter().zip(&base.pressure_pa) {
            assert!((pressure - offset - original).abs() < 2e-10);
        }
        let (_, rest) = channel(4, 0, offset, offset);
        let rest = simple_flow(&domain, &fluid(), &rest, &gate).unwrap();
        assert!(rest.velocity_m_s.iter().all(|v| *v == [0.0; 3]));
        assert!(rest.pressure_pa.iter().all(|p| *p == offset));
        assert_eq!(rest.report.inflow_m3_s, 0.0);
        assert_eq!(rest.report.outflow_m3_s, 0.0);
    }
    let mut zero_outlet = config;
    zero_outlet.faces[1] = FvBoundary::Outlet;
    let flow = simple_flow(&domain, &fluid(), &zero_outlet, &gate).unwrap();
    assert_eq!(flow, base);
}

#[test]
fn pressure_opening_flow_reverses_on_every_axis_and_respects_material_viscosity() {
    let gate = CancelGate::new();
    let mut throughputs = Vec::new();
    for axis in 0..3 {
        for (left, right, sign) in [(0.5, 0.0, 1.0), (0.0, 0.5, -1.0)] {
            let (domain, config) = channel(4, axis, left, right);
            let mut rates = Vec::new();
            for viscosity in [0.25, 0.5] {
                let flow = simple_flow_with_viscosity(
                    &domain,
                    &fluid(),
                    &vec![viscosity; domain.cell_count()],
                    &config,
                    &gate,
                )
                .unwrap();
                let q = flow
                    .field
                    .boundary_outflow(&domain, Face3::ALL[2 * axis + 1]);
                assert!(q * sign > 0.0);
                assert!(flow.report.max_divergence_m3_s < 1e-12);
                assert!((flow.report.inflow_m3_s - q.abs()).abs() < 1e-12);
                rates.push(q.abs());
            }
            assert!(
                (rates[0] / rates[1] - 2.0).abs() < 2e-7,
                "axis {axis}: {rates:?}"
            );
            throughputs.push(rates[0]);
        }
    }
    assert!(
        throughputs
            .iter()
            .all(|q| (q / throughputs[0] - 1.0).abs() < 2e-7)
    );
}

#[test]
fn pressure_opening_backflow_uses_the_existing_thermal_reservoir() {
    let gate = CancelGate::new();
    for (left, right, temperature) in [(0.5, 0.0, 320.0), (0.0, 0.5, 280.0)] {
        let (domain, config) = channel(4, 0, left, right);
        let flow = simple_flow(&domain, &fluid(), &config, &gate).unwrap();
        let mut faces = [ThermalFace::Adiabatic; 6];
        faces[0] = ThermalFace::Outflow {
            backflow_temperature: 320.0,
        };
        faces[1] = ThermalFace::Outflow {
            backflow_temperature: 280.0,
        };
        let energy = solve_energy(
            &domain,
            &fluid(),
            &[],
            &flow.field,
            &ThermalSetup::new(faces),
            &EnergyConfig::default(),
            &gate,
        )
        .unwrap();
        assert!(
            energy
                .temperature
                .iter()
                .all(|t| (t - temperature).abs() < 1e-6)
        );
        assert!(energy.report.balance.residual_w.abs() < 1e-8);
    }
}

#[test]
fn pressure_opening_drives_the_original_transient_momentum_owner() {
    let gate = CancelGate::new();
    let (domain, config) = channel(4, 0, 0.5, 0.0);
    let steady = simple_flow(&domain, &fluid(), &config, &gate).unwrap();
    let mut time = UnsteadyConfig::new(0.1, 12);
    time.inner_tolerance = 1e-9;
    time.inner_iterations = 300;
    let started = simple_unsteady(&domain, &fluid(), &config, &time, |_| 1.0, &gate).unwrap();
    let energies: Vec<_> = started.records.iter().map(|r| r.kinetic_energy_j).collect();
    assert!(energies[0] > 0.0);
    assert!(energies.windows(2).all(|pair| pair[1] > pair[0]));
    let q = started.flow.field.boundary_outflow(&domain, Face3::XMax);
    assert!(q > 0.9 * steady.report.outflow_m3_s && q < steady.report.outflow_m3_s);
    assert!(started.flow.report.max_divergence_m3_s < 1e-12);
    let (_, shifted) = channel(4, 0, 1024.5, 1024.0);
    let shifted = simple_unsteady(&domain, &fluid(), &shifted, &time, |_| 1.0, &gate).unwrap();
    assert_eq!(started.records, shifted.records);
    assert_eq!(started.flow.velocity_m_s, shifted.flow.velocity_m_s);
}

fn pressure_plug(axis: usize) -> (VoxelDomain, SimpleConfig) {
    let (domain, mut config) = channel(1, axis, 0.5, 0.0);
    let transverse = (axis + 1) % 3;
    config.faces[2 * transverse] = FvBoundary::Symmetry;
    config.faces[2 * transverse + 1] = FvBoundary::Symmetry;
    (domain, config)
}

#[test]
fn scheduled_pressure_retains_native_be_bdf2_inertia_and_reverses_on_all_axes() {
    let gate = CancelGate::new();
    for axis in 0..3 {
        for scheme in [TimeScheme::BackwardEuler, TimeScheme::Bdf2] {
            let (domain, config) = pressure_plug(axis);
            let mut time = UnsteadyConfig::new(0.125, 4);
            time.scheme = scheme;
            time.inner_tolerance = 1e-10;
            time.inner_iterations = 400;
            time.probe = Some(0);
            let flow = simple_unsteady_with_pressure_schedule(
                &domain,
                &fluid(),
                &config,
                &time,
                |_| 1.0,
                |t| if t <= 0.25 { 1.0 } else { -4.0 },
                &gate,
            )
            .unwrap();
            // Symmetry sides admit an exact uniform plug. Its native pressure
            // gradient supplies du/dt = dp/(rho L); convection and shear vanish.
            // This recurrence checks the actual solved face field and inertia.
            let mut previous = 0.0;
            let mut older = 0.0;
            for (step, record) in flow.records.iter().enumerate() {
                let factor = if step < 2 { 1.0 } else { -4.0 };
                let impulse = time.time_step_s * (0.5 / 2.0) * factor;
                let expected = match scheme {
                    TimeScheme::Bdf2 if step > 0 => (2.0 * previous - 0.5 * older + impulse) / 1.5,
                    _ => previous + impulse,
                };
                let mut velocity = [0.0; 3];
                velocity[axis] = expected;
                for (actual, expected) in record.probe_velocity_m_s.iter().zip(velocity) {
                    assert!(
                        (actual - expected).abs() < 2e-9,
                        "axis={axis}, scheme={scheme:?}, step={step}: {record:?}"
                    );
                }
                // rho=1, physical volume=2: K = rho V u^2/2 = u^2,
                // including the original opening half-cell masses.
                assert!((record.kinetic_energy_j - expected * expected).abs() < 2e-9);
                older = previous;
                previous = expected;
            }
            assert!(previous < 0.0);
            let q = flow
                .flow
                .field
                .boundary_outflow(&domain, Face3::ALL[2 * axis + 1]);
            assert!((q - previous).abs() < 2e-9);
            assert!(flow.flow.report.max_divergence_m3_s < 1e-10);
            for c in 0..domain.cell_count() {
                let at = domain.coords(c);
                let x = domain.center(at[0], at[1], at[2])[axis];
                let pressure = 0.5 - (-4.0) * (0.5 / 2.0) * x;
                assert!((flow.flow.pressure_pa[c] - pressure).abs() < 2e-8);
            }
        }
    }
}

#[test]
fn oscillatory_pressure_native_flow_refines_in_time_and_conserves_volume() {
    let gate = CancelGate::new();
    let (domain, config) = pressure_plug(0);
    let omega = core::f64::consts::TAU;
    for scheme in [TimeScheme::BackwardEuler, TimeScheme::Bdf2] {
        let mut errors = Vec::new();
        for steps in [24, 48, 96] {
            let mut time = UnsteadyConfig::new(0.75 / steps as f64, steps);
            time.scheme = scheme;
            time.inner_tolerance = 1e-11;
            time.inner_iterations = 400;
            time.probe = Some(0);
            let flow = simple_unsteady_with_pressure_schedule(
                &domain,
                &fluid(),
                &config,
                &time,
                |_| 1.0,
                |t| (omega * t).cos(),
                &gate,
            )
            .unwrap();
            // du/dt = (dp/rho L) cos(omega t), u(0)=0. At t=3/4
            // period, u=-dp/(rho L omega), after actual flow reversal.
            let expected = -0.25 / omega;
            let actual = flow.records.last().unwrap().probe_velocity_m_s[0];
            errors.push((actual - expected).abs());
            assert!(flow.records.iter().any(|r| r.probe_velocity_m_s[0] > 0.01));
            assert!(actual < 0.0);
            assert!(flow.flow.report.max_divergence_m3_s < 1e-10);
            assert!(flow.mean_field.max_divergence(&domain) < 1e-10);
        }
        let minimum_ratio = match scheme {
            TimeScheme::BackwardEuler => 1.9,
            TimeScheme::Bdf2 => 3.3,
        };
        assert!(
            errors.windows(2).all(|e| e[0] / e[1] > minimum_ratio),
            "scheme={scheme:?}: {errors:?}"
        );
    }
}

#[test]
fn pressure_schedule_keeps_profiles_component_gauges_and_constant_schedule_parity() {
    let gate = CancelGate::new();
    let (domain, config) = separated_pressure_channels(0, [0.0; 2], [0.5, -0.25]);
    let (_, shifted) =
        separated_pressure_channels(0, [1_048_576.0, -1_099_511_627_776.0], [0.5, -0.25]);
    let mut time = UnsteadyConfig::new(0.1, 4);
    time.inner_tolerance = 1e-9;
    time.inner_iterations = 400;
    let original = simple_unsteady(&domain, &fluid(), &config, &time, |_| 1.0, &gate).unwrap();
    let constant = simple_unsteady_with_pressure_schedule(
        &domain,
        &fluid(),
        &config,
        &time,
        |_| 1.0,
        |_| 1.0,
        &gate,
    )
    .unwrap();
    assert_eq!(constant, original);
    let samples = core::cell::RefCell::new(Vec::new());
    let schedule = |t| if t <= 0.2 { 1.0 } else { -2.0 };
    let driven = simple_unsteady_with_pressure_schedule(
        &domain,
        &fluid(),
        &config,
        &time,
        |_| 1.0,
        |t| {
            samples.borrow_mut().push(t);
            schedule(t)
        },
        &gate,
    )
    .unwrap();
    let offset = simple_unsteady_with_pressure_schedule(
        &domain,
        &fluid(),
        &shifted,
        &time,
        |_| 1.0,
        schedule,
        &gate,
    )
    .unwrap();
    assert_eq!(
        *samples.borrow(),
        (1..=time.steps)
            .map(|step| step as f64 * time.time_step_s)
            .collect::<Vec<_>>()
    );
    assert_eq!(driven.records, offset.records);
    assert_eq!(driven.flow.velocity_m_s, offset.flow.velocity_m_s);
    assert_eq!(driven.flow.field, offset.flow.field);
    assert_eq!(driven.flow.report, offset.flow.report);
    assert_eq!(driven.mean_field, offset.mean_field);
    assert!(driven.flow.report.max_divergence_m3_s < 1e-10);
}

#[test]
fn pressure_schedule_refuses_nonfinite_overflow_and_cancellation_without_mutating_inputs() {
    let gate = CancelGate::new();
    let (domain, config) = channel(1, 0, 4.0, 0.0);
    let original = config.clone();
    let mut time = UnsteadyConfig::new(0.1, 2);
    time.inner_tolerance = 1e-9;
    time.inner_iterations = 400;
    for (sample, field) in [
        (f64::NAN, "unsteady.pressure_schedule"),
        (f64::INFINITY, "unsteady.pressure_schedule"),
        (f64::NEG_INFINITY, "unsteady.pressure_schedule"),
        (f64::MAX, "unsteady.pressure_difference_pa"),
    ] {
        let result = simple_unsteady_with_pressure_schedule(
            &domain,
            &fluid(),
            &config,
            &time,
            |_| 1.0,
            |t| if t <= 0.1 { 1.0 } else { sample },
            &gate,
        );
        assert!(matches!(result, Err(ChtError::InvalidInput { field: got, .. }) if got == field));
        assert_eq!(config, original);
    }
    let result = simple_unsteady_with_pressure_schedule(
        &domain,
        &fluid(),
        &config,
        &time,
        |_| 1.0,
        |_| {
            gate.request();
            1.0
        },
        &gate,
    );
    assert!(matches!(result, Err(ChtError::Cancelled)));
    assert_eq!(config, original);
}

#[test]
fn pressure_opening_refuses_nonfinite_and_unrepresentable_boundary_differences() {
    let gate = CancelGate::new();
    for pressure in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let (domain, config) = channel(4, 0, pressure, 0.0);
        assert!(matches!(
            simple_flow(&domain, &fluid(), &config, &gate),
            Err(ChtError::InvalidInput {
                field: "simple.pressure_pa",
                ..
            })
        ));
    }
    let (domain, config) = channel(4, 0, f64::MAX, -f64::MAX);
    assert!(matches!(
        simple_flow(&domain, &fluid(), &config, &gate),
        Err(ChtError::InvalidInput {
            field: "simple.pressure_difference_pa",
            ..
        })
    ));
    let (_, config) = channel(4, 0, 0.5, 0.0);
    gate.request();
    assert!(matches!(
        simple_flow(&domain, &fluid(), &config, &gate),
        Err(ChtError::Cancelled)
    ));
}

// One actual MAC momentum predictor, with a resolved normal carrier frozen
// while solving its transverse component. Independent three-row equations
// below include the declared relaxation; this is not a pressure-flow mock.
fn opening_transverse_predictor(axis: usize, reverse: bool, counterflow: bool) -> [f64; 3] {
    let tangent = (axis + 1) % 3;
    let mut dims = [1; 3];
    dims[axis] = 3;
    dims[tangent] = 2;
    let domain = VoxelDomain::new(dims[0], dims[1], dims[2], 1.0).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[2 * axis] = FvBoundary::Pressure { pressure_pa: 0.0 };
    faces[2 * axis + 1] = FvBoundary::Pressure { pressure_pa: 0.0 };
    let source = 2 * axis + usize::from(reverse);
    // Outlet is the same pressure opening, including when it receives flow.
    faces[source] = FvBoundary::Outlet;
    let mut config = SimpleConfig::new(faces);
    config.scheme = ConvectionScheme::Upwind;
    config.velocity_relaxation = 0.5;
    config.momentum_tolerance = 1e-12;
    let mut external = [0.0; 3];
    external[tangent] = 0.25;
    config.backflow_tangential_velocity_m_s[source] = Some(external);
    external[tangent] = -7.0;
    config.backflow_tangential_velocity_m_s[source ^ 1] = Some(external);
    let fluid = fluid();
    admit(&domain, &fluid, &config).unwrap();
    let mut solver = Solver::new(&domain, &fluid, &config);
    let sign = if reverse { -1.0 } else { 1.0 };
    if counterflow {
        // Equal incoming/outgoing half-faces: net flux is zero, but incoming
        // reservoir momentum and outgoing interior momentum are both real.
        for row in 0..2 {
            let mut face = [0; 3];
            face[axis] = if reverse { 3 } else { 0 };
            face[tangent] = row;
            let index = solver.comps[axis].index(face);
            solver.vel[axis][index] = if row == 0 { sign } else { -sign };
        }
    } else {
        solver.vel[axis].fill(sign);
    }
    solver.momentum(tangent, true, &CancelGate::new()).unwrap();
    std::array::from_fn(|position| {
        let mut face = [0; 3];
        face[axis] = if reverse { 2 - position } else { position };
        face[tangent] = 1;
        solver.vel[tangent][solver.comps[tangent].index(face)]
    })
}

#[test]
fn pressure_opening_backflow_momentum_obeys_each_half_face_budget_on_all_axes() {
    for axis in 0..3 {
        for reverse in [false, true] {
            for counterflow in [false, true] {
                // rho=1, mu=1/4, unit geometry, alpha=1/2. Ordinary inflow
                // carries 1/4 N; opposed half-flows carry 1/8 N even at net Q=0.
                let (diagonal, lower, upper, incoming, outgoing_rate) = if counterflow {
                    ([2.5, 2.0, 1.5], 0.25, 0.25, 0.125, 0.5)
                } else {
                    ([3.5, 4.0, 3.5], 1.25, 0.25, 0.25, 1.0)
                };
                let middle = (lower * incoming / diagonal[0])
                    / (diagonal[1] - lower * upper / diagonal[0] - upper * lower / diagonal[2]);
                let expected = [
                    (incoming + upper * middle) / diagonal[0],
                    middle,
                    lower * middle / diagonal[2],
                ];
                let actual = opening_transverse_predictor(axis, reverse, counterflow);
                for (&a, &b) in actual.iter().zip(&expected) {
                    assert!(
                        (a - b).abs() < 1e-12,
                        "axis={axis}, reverse={reverse}, counterflow={counterflow}: {actual:?} != {expected:?}"
                    );
                }
                // Sum the independently represented boundary momentum, normal
                // viscous resistance and numerical under-relaxation terms.
                let outgoing = outgoing_rate * actual[if counterflow { 0 } else { 2 }];
                let viscous = 0.5 * actual.iter().sum::<f64>();
                let relaxation: f64 = actual
                    .iter()
                    .zip(diagonal)
                    .map(|(&velocity, diagonal)| 0.5 * diagonal * velocity)
                    .sum();
                assert!((incoming - outgoing - viscous - relaxation).abs() < 1e-12);
            }
        }
    }
}

#[test]
fn pressure_opening_backflow_momentum_reaches_steady_and_transient_native_flow() {
    let gate = CancelGate::new();
    let (domain, mut config) = channel(4, 0, 0.5, 0.0);
    let baseline = simple_flow(&domain, &fluid(), &config, &gate).unwrap();
    config.backflow_tangential_velocity_m_s[0] = Some([0.0; 3]);
    config.backflow_tangential_velocity_m_s[1] = Some([0.0; 3]);
    assert_eq!(
        simple_flow(&domain, &fluid(), &config, &gate).unwrap(),
        baseline
    );
    config.backflow_tangential_velocity_m_s[0] = Some([0.0, 0.03, 0.0]);
    let driven = simple_flow(&domain, &fluid(), &config, &gate).unwrap();
    assert!(driven.velocity_m_s.iter().any(|v| v[1].abs() > 1e-6));
    assert!(driven.report.max_divergence_m3_s < 1e-10);
    let mut time = UnsteadyConfig::new(0.1, 3);
    time.inner_tolerance = 1e-8;
    time.inner_iterations = 300;
    let started = simple_unsteady(&domain, &fluid(), &config, &time, |_| 1.0, &gate).unwrap();
    assert!(started.flow.velocity_m_s.iter().any(|v| v[1].abs() > 1e-7));
    assert!(started.flow.report.max_divergence_m3_s < 1e-10);
}

#[test]
fn pressure_opening_backflow_momentum_refuses_unsupported_components_and_work_model() {
    let (domain, original) = channel(4, 0, 0.5, 0.0);
    for (face, velocity) in [
        (0, [0.1, 0.0, 0.0]),
        (0, [0.0, f64::NAN, 0.0]),
        (0, [0.0, f64::INFINITY, 0.0]),
        (2, [0.1, 0.0, 0.0]),
    ] {
        let mut config = original.clone();
        config.backflow_tangential_velocity_m_s[face] = Some(velocity);
        assert!(matches!(
            admit(&domain, &fluid(), &config),
            Err(ChtError::InvalidInput {
                field: "simple.backflow_tangential_velocity_m_s",
                ..
            })
        ));
    }
    let mut prescribed = original.clone();
    prescribed.faces[0] = FvBoundary::Inlet {
        velocity: [0.1, 0.0, 0.0],
    };
    prescribed.backflow_tangential_velocity_m_s[0] = Some([0.0, 0.1, 0.0]);
    assert!(admit(&domain, &fluid(), &prescribed).is_err());
    let controls = viscous_heat::ViscousHeatControls {
        maximum_reynolds: 0.1,
        power_tolerance_w: 1e-8,
    };
    let mut config = original;
    config.backflow_tangential_velocity_m_s[0] = Some([0.0; 3]);
    assert!(
        viscous_heat::pressure::admit_pressure_work_heating(&domain, &config, controls).is_ok()
    );
    config.backflow_tangential_velocity_m_s[0] = Some([0.0, 0.03, 0.0]);
    assert!(
        viscous_heat::pressure::admit_pressure_work_heating(&domain, &config, controls).is_err()
    );
}

#[test]
fn pressure_projection_polishes_quiet_component_before_caloric_transport() {
    use crate::conjugate::{
        AdvectiveEnthalpyField, AffineSpecificEnthalpy, CellSink, SolidMaterial, Voxel,
    };

    // Five fluid cells: one drained cell and a separated 2x2 recirculating
    // pocket. A single global pressure solve can accept a large LOCAL defect
    // in that pocket because its flux is tiny relative to the other component.
    let mut domain = VoxelDomain::new(4, 2, 1, 1.0).unwrap();
    for [x, y] in [[1, 0], [1, 1], [0, 1]] {
        domain.set(x, y, 0, Voxel::Solid(0));
    }
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    let mut config = SimpleConfig::new(faces);
    config.pressure_solver = PressureSolver::IluBicgstab;
    let gate = CancelGate::new();
    let fluid = fluid();
    let epsilon = 1e-12;
    let seed = |solver: &mut Solver<'_>| {
        for axis in 0..3 {
            for &face in &solver.comps[axis].unknowns {
                solver.d[axis][face] = 1.0;
            }
        }
        for (axis, at, speed) in [
            (0, [0, 0, 0], 1.0),
            (0, [3, 0, 0], 2.0 * epsilon),
            (0, [3, 1, 0], -epsilon),
            (1, [2, 1, 0], -epsilon),
            (1, [3, 1, 0], epsilon),
        ] {
            let face = solver.comps[axis].index(at);
            solver.vel[axis][face] = speed;
        }
    };
    let mut single = Solver::new(&domain, &fluid, &config);
    seed(&mut single);
    single.correct(1e-8, &gate).unwrap();
    // The actual ILU(0) first correction is epsilon*[-5/8,1/8,-3/8,-1/16]
    // in the pocket. Its global true residual is ~1.13e-13, but the worst
    // cell retains imbalance 3*epsilon/16 and relative continuity defect 3/43.
    let field = single.field();
    let cell = domain.index(3, 0, 0);
    let magnitude: f64 = (0..6)
        .map(|face| field.outward(&domain, cell, face).abs())
        .sum();
    assert!((field.net_outflow(&domain, cell).abs() / magnitude - 3.0 / 43.0).abs() < 1e-12);

    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.faces[0] = ThermalFace::Outflow {
        backflow_temperature: 300.0,
    };
    thermal.faces[1] = ThermalFace::Outflow {
        backflow_temperature: 300.0,
    };
    thermal.cell_sinks.push(CellSink {
        cell: 0,
        conductance_w_k: 1.0,
        temperature_k: 300.0,
    });
    let caloric: Vec<_> = (0..domain.cell_count())
        .map(|cell| {
            domain.is_fluid(cell).then_some(AffineSpecificEnthalpy {
                temperature_k: 300.0,
                enthalpy_j_kg: 0.0,
                heat_capacity_j_kg_k: 2.0,
            })
        })
        .collect();
    let heat = |field: &FlowField| {
        EnergyConfig {
            scheme: ConvectionScheme::Upwind,
            ..EnergyConfig::default()
        }
        .solve_with_enthalpy(
            &domain,
            &fluid,
            &[SolidMaterial::new("partition", 1.0)],
            field,
            &thermal,
            &vec![[1.0; 3]; domain.cell_count()],
            AdvectiveEnthalpyField {
                cells: &caloric,
                boundary_j_kg: [Some(0.0), Some(0.0), None, None, None, None],
            },
            &gate,
        )
    };
    assert!(matches!(
        heat(&field),
        Err(ChtError::InvalidInput {
            field: "thermal.advective_enthalpy",
            ..
        })
    ));

    for pressure_solver in [PressureSolver::IluBicgstab, PressureSolver::AmgCg] {
        config.pressure_solver = pressure_solver;
        let mut solver = Solver::new(&domain, &fluid, &config);
        seed(&mut solver);
        solver.project(&gate).unwrap();
        let velocities = solver.vel.clone();
        let pressure = solver.pressure.clone();
        let iterations = solver.pressure_krylov;
        solver.project(&gate).unwrap();
        assert_eq!(solver.vel, velocities);
        assert_eq!(solver.pressure, pressure);
        assert_eq!(solver.pressure_krylov, iterations);
        let field = solver.field();
        for cell in 0..domain.cell_count() {
            let magnitude: f64 = (0..6)
                .map(|face| field.outward(&domain, cell, face).abs())
                .sum();
            assert!(field.net_outflow(&domain, cell).abs() <= 128.0 * f64::EPSILON * magnitude);
        }
        // Pressure corrects the native faces; it does not clip the quiet
        // circulation to rest to manufacture continuity.
        assert!(field.flux_x(3, 0, 0).abs() > epsilon);
        let energy = heat(&field).unwrap();
        assert!(
            energy
                .temperature
                .iter()
                .all(|temperature| (*temperature - 300.0).abs() < 1e-8)
        );
        assert!(energy.report.balance.residual_w.abs() < 1e-8);
    }
}

#[test]
fn pressure_projection_refuses_incompatible_fixed_flux_without_clipping_it() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Inlet {
        velocity: [1.0, 0.0, 0.0],
    };
    let config = SimpleConfig::new(faces);
    let mut solver = Solver::new(&domain, &fluid(), &config);
    let original = solver.vel.clone();
    assert!(matches!(
        solver.project(&CancelGate::new()),
        Err(ChtError::SolverNotConverged {
            system: "pressure continuity",
            ..
        })
    ));
    assert_eq!(solver.vel, original);
    let pressure = solver.pressure.clone();
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(solver.project(&gate), Err(ChtError::Cancelled)));
    assert_eq!(solver.vel, original);
    assert_eq!(solver.pressure, pressure);
}

#[cfg(feature = "enthalpy-coupling")]
#[test]
fn pressure_opening_native_momentum_restart_retains_its_numerical_datum() {
    let (domain, config) = channel(4, 0, 1024.5, 1024.0);
    let fluid = fluid();
    let gate = CancelGate::new();
    let mut solver = Solver::new(&domain, &fluid, &config);
    for _ in 0..5 {
        solver.sweep(&gate).unwrap();
    }
    let state = solver.momentum_state();
    let mut restored = Solver::new(&domain, &fluid, &config);
    restored.restore_momentum(&state).unwrap();
    assert_eq!(restored.pressure, solver.pressure);
    assert_eq!(restored.pressure_reference_pa, solver.pressure_reference_pa);
    assert_eq!(restored.vel, solver.vel);
}

fn plane_profile(
    domain: &VoxelDomain,
    side: usize,
    pressure: impl Fn([f64; 3]) -> f64,
) -> Vec<f64> {
    let axis = side / 2;
    let n = domain.dims();
    let plane = if side % 2 == 0 { 0 } else { n[axis] - 1 };
    let boundary = if side % 2 == 0 {
        0.0
    } else {
        (0..n[axis]).map(|i| domain.width(axis, i)).sum()
    };
    (0..domain.cell_count())
        .filter_map(|cell| {
            let at = domain.coords(cell);
            (at[axis] == plane).then(|| {
                let mut point = domain.center(at[0], at[1], at[2]);
                point[axis] = boundary;
                pressure(point)
            })
        })
        .collect()
}

#[test]
fn spatial_pressure_uniform_profiles_retain_native_results_on_every_axis() {
    let gate = CancelGate::new();
    for axis in 0..3 {
        for (inlet, outlet) in [(0.5, 0.0), (0.0, 0.5)] {
            let (domain, mut config) = channel(4, axis, inlet, outlet);
            let uniform = simple_flow(&domain, &fluid(), &config, &gate).unwrap();
            config.pressure_profiles_pa[2 * axis] =
                Some(plane_profile(&domain, 2 * axis, |_| inlet));
            config.pressure_profiles_pa[2 * axis + 1] =
                Some(plane_profile(&domain, 2 * axis + 1, |_| outlet));
            // Profiles replace the scalar pressure, including the datum.
            config.faces[2 * axis] = FvBoundary::Pressure { pressure_pa: -99.0 };
            config.faces[2 * axis + 1] = FvBoundary::Outlet;
            let resolved = simple_flow(&domain, &fluid(), &config, &gate).unwrap();
            assert_eq!(resolved, uniform);
        }
    }
}

#[test]
fn spatial_pressure_balances_hydrostatic_force_without_changing_channel_flow() {
    let gate = CancelGate::new();
    let fluid = fluid();
    for axis in 0..3 {
        let gravity_axis = (axis + 1) % 3;
        let wall_axis = (axis + 2) % 3;
        let mut widths = [vec![1.0], vec![1.0], vec![1.0]];
        widths[axis] = vec![0.25, 0.5, 0.75, 0.5];
        widths[gravity_axis] = vec![0.25, 0.75, 1.0];
        widths[wall_axis] = vec![0.25; 4];
        let domain = VoxelDomain::graded(widths).unwrap();
        let mut faces = [FvBoundary::Symmetry; 6];
        faces[2 * axis] = FvBoundary::Pressure { pressure_pa: 0.5 };
        faces[2 * axis + 1] = FvBoundary::Outlet;
        faces[2 * wall_axis] = FvBoundary::wall();
        faces[2 * wall_axis + 1] = FvBoundary::wall();
        let mut config = SimpleConfig::new(faces);
        config.tolerance = 1e-10;
        config.pressure_tolerance = 1e-12;
        let mu = vec![0.25; domain.cell_count()];
        let uniform = simple_flow_with_viscosity(&domain, &fluid, &mu, &config, &gate).unwrap();
        let pressure = |point: [f64; 3]| 16.0 - 0.25 * point[axis] - 0.25 * point[gravity_axis];
        for side in [2 * axis, 2 * axis + 1] {
            config.pressure_profiles_pa[side] = Some(plane_profile(&domain, side, pressure));
        }
        let mut force = [0.0; 3];
        force[gravity_axis] = -0.25;
        let force = vec![force; domain.cell_count()];
        let loaded =
            simple_flow_with_viscosity_and_force(&domain, &fluid, &mu, &force, &config, &gate)
                .unwrap();
        for (actual, reference) in loaded.velocity_m_s.iter().zip(&uniform.velocity_m_s) {
            for component in 0..3 {
                assert!((actual[component] - reference[component]).abs() < 2e-8);
            }
        }
        for (cell, &actual) in loaded.pressure_pa.iter().enumerate() {
            let at = domain.coords(cell);
            assert!((actual - pressure(domain.center(at[0], at[1], at[2]))).abs() < 2e-7);
        }
        // Pure hydrostatics is an exact fixed point of the same discrete
        // momentum equations. Seed the analytic state to isolate the spatial
        // force balance from a near-zero-flow relative stopping criterion.
        let hydrostatic = |point: [f64; 3]| 16.0 - 0.25 * point[gravity_axis];
        for side in [2 * axis, 2 * axis + 1] {
            config.pressure_profiles_pa[side] = Some(plane_profile(&domain, side, hydrostatic));
        }
        admit(&domain, &fluid, &config).unwrap();
        let mut rest = Solver::new(&domain, &fluid, &config);
        rest.set_material_transport(&mu, Some(&force), &gate)
            .unwrap();
        for cell in 0..domain.cell_count() {
            let at = domain.coords(cell);
            rest.pressure[cell] =
                hydrostatic(domain.center(at[0], at[1], at[2])) - rest.pressure_reference_pa.at(cell);
        }
        for component in 0..3 {
            assert_eq!(rest.momentum(component, false, &gate).unwrap().rms, 0.0);
        }
        assert!(rest.vel.iter().flatten().all(|&velocity| velocity == 0.0));
    }
}

#[test]
fn spatial_pressure_drives_opposed_channel_streams_and_ignores_solid_faces() {
    use crate::conjugate::Voxel;

    let mut domain = VoxelDomain::graded([vec![0.25; 8], vec![0.5; 5], vec![1.0]]).unwrap();
    for x in 0..8 {
        domain.set(x, 2, 0, Voxel::Solid(0));
    }
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    faces[2] = FvBoundary::wall();
    faces[3] = FvBoundary::wall();
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-10;
    config.pressure_tolerance = 1e-12;
    config.pressure_profiles_pa[0] = Some(vec![0.5, 0.5, 1e12, -0.25, -0.25]);
    config.pressure_profiles_pa[1] = Some(vec![0.0, 0.0, -1e12, 0.0, 0.0]);
    let solved = simple_flow(&domain, &fluid(), &config, &CancelGate::new()).unwrap();
    // Each channel has H=1, L=2, mu=1/4 and two transverse cells:
    // Q=dp H^3/(12 mu L)*(1+2/2^2), with each pressure sign retained.
    let low = solved.field.flux_x(8, 0, 0) + solved.field.flux_x(8, 1, 0);
    let high = solved.field.flux_x(8, 3, 0) + solved.field.flux_x(8, 4, 0);
    assert!((low - 0.125).abs() < 1e-8, "{low}");
    assert!((high + 0.0625).abs() < 1e-8, "{high}");
    assert_eq!(solved.field.flux_x(0, 2, 0), 0.0);
    assert_eq!(solved.field.flux_x(8, 2, 0), 0.0);
}

#[test]
fn spatial_pressure_leading_solid_face_cannot_erase_the_active_driving_head() {
    use crate::conjugate::Voxel;

    let mut domain =
        VoxelDomain::graded([vec![0.25; 8], vec![0.5; 3], vec![1.0]]).unwrap();
    for x in 0..8 {
        domain.set(x, 0, 0, Voxel::Solid(0));
    }
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Pressure { pressure_pa: 1.0 };
    faces[1] = FvBoundary::Outlet;
    faces[2] = FvBoundary::wall();
    faces[3] = FvBoundary::wall();
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-10;
    config.pressure_tolerance = 1e-12;
    let gate = CancelGate::new();
    let reference = simple_flow(&domain, &fluid(), &config, &gate).unwrap();
    assert!((reference.report.outflow_m3_s - 0.25).abs() < 1e-8);
    for unused in [1e18, f64::MAX] {
        config.pressure_profiles_pa[0] = Some(vec![unused, 1.0, 1.0]);
        config.pressure_profiles_pa[1] = Some(vec![-unused, 0.0, 0.0]);
        // A 1 Pa active head must survive arbitrary FINITE inactive values.
        // f64::MAX additionally proves unused doubled pressure differences
        // are never admitted as if they loaded a fluid face.
        let actual = simple_flow(&domain, &fluid(), &config, &gate).unwrap();
        assert_eq!(actual, reference);
    }
    config.pressure_profiles_pa[0] = Some(vec![f64::NAN, 1.0, 1.0]);
    assert!(matches!(
        simple_flow(&domain, &fluid(), &config, &gate),
        Err(ChtError::InvalidInput { field: "simple.pressure_profiles_pa", .. })
    ));
}

#[test]
fn spatial_pressure_work_uses_each_original_face_and_closes_physical_heat() {
    use super::viscous_heat::{ViscousHeatControls, creeping_flow_with_pressure_work_heat};

    let domain = VoxelDomain::graded([vec![0.5, 1.0, 1.5], vec![0.75, 1.25], vec![1.0]]).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    faces[2] = FvBoundary::wall();
    faces[3] = FvBoundary::wall();
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-12;
    config.pressure_tolerance = 1e-12;
    config.momentum_tolerance = 1e-10;
    config.pressure_profiles_pa[0] = Some(vec![2.0, 6.0]);
    config.pressure_profiles_pa[1] = Some(vec![0.0; 2]);
    let controls = ViscousHeatControls {
        maximum_reynolds: 0.1,
        power_tolerance_w: 1e-7,
    };
    let mu = vec![100.0; domain.cell_count()];
    let gate = CancelGate::new();
    let (flow, heat) =
        creeping_flow_with_pressure_work_heat(&domain, &fluid(), &mu, &config, controls, &gate)
            .unwrap();
    let supply = 2.0 * flow.field.flux_x(0, 0, 0) + 6.0 * flow.field.flux_x(0, 1, 0);
    assert!(supply > 0.0);
    assert!((supply - heat.supplied_pressure_power_w).abs() < 1e-10);
    assert!((supply - heat.heat_w).abs() < controls.power_tolerance_w);
    // The area-mean pressure times the integrated flow loses the covariance.
    let mean = (0.75 * 2.0 + 1.25 * 6.0) / 2.0;
    let averaged = mean * (flow.field.flux_x(0, 0, 0) + flow.field.flux_x(0, 1, 0));
    assert!((averaged - supply).abs() > 1e-5 * supply);
    for profile in config.pressure_profiles_pa.iter_mut().flatten() {
        for pressure in profile {
            *pressure += 1024.0;
        }
    }
    let (shifted, shifted_heat) =
        creeping_flow_with_pressure_work_heat(&domain, &fluid(), &mu, &config, controls, &gate)
            .unwrap();
    assert_eq!(shifted.velocity_m_s, flow.velocity_m_s);
    assert_eq!(shifted_heat.cell_heat_w, heat.cell_heat_w);
    assert_eq!(shifted_heat.pressure_power_w, heat.pressure_power_w);
    assert!(shifted_heat.residual_w.abs() < controls.power_tolerance_w);
}

#[test]
fn spatial_pressure_refuses_incompatible_shape_role_values_and_cancellation() {
    let (domain, config) = channel(4, 0, 0.5, 0.0);
    let gate = CancelGate::new();
    for values in [vec![], vec![0.0; 3], vec![0.0, 0.0, f64::NAN, 0.0]] {
        let mut invalid = config.clone();
        invalid.pressure_profiles_pa[0] = Some(values);
        assert!(matches!(
            simple_flow(&domain, &fluid(), &invalid, &gate),
            Err(ChtError::InvalidInput {
                field: "simple.pressure_profiles_pa",
                ..
            })
        ));
    }
    let mut invalid = config.clone();
    invalid.pressure_profiles_pa[2] = Some(vec![0.0; 8]);
    assert!(matches!(
        simple_flow(&domain, &fluid(), &invalid, &gate),
        Err(ChtError::InvalidInput {
            field: "simple.pressure_profiles_pa",
            ..
        })
    ));
    invalid = config.clone();
    invalid.pressure_profiles_pa[0] = Some(vec![-f64::MAX, f64::MAX, 0.0, 0.0]);
    assert!(matches!(
        simple_flow(&domain, &fluid(), &invalid, &gate),
        Err(ChtError::InvalidInput {
            field: "simple.pressure_difference_pa",
            ..
        })
    ));
    let mut valid = config;
    valid.pressure_profiles_pa[0] = Some(vec![0.5; 4]);
    gate.request();
    assert!(matches!(
        simple_flow(&domain, &fluid(), &valid, &gate),
        Err(ChtError::Cancelled)
    ));
}


fn separated_pressure_channels(
    axis: usize,
    offsets: [f64; 2],
    heads: [f64; 2],
) -> (VoxelDomain, SimpleConfig) {
    use crate::conjugate::Voxel;

    let split = (axis + 1) % 3;
    let mut widths = [vec![1.0], vec![1.0], vec![1.0]];
    widths[axis] = vec![0.25; 8];
    widths[split] = vec![0.5; 5];
    let mut domain = VoxelDomain::graded(widths).unwrap();
    for cell in 0..domain.cell_count() {
        let at = domain.coords(cell);
        if at[split] == 2 {
            domain.set(at[0], at[1], at[2], Voxel::Solid(0));
        }
    }
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[2 * axis] = FvBoundary::Outlet;
    faces[2 * axis + 1] = FvBoundary::Outlet;
    faces[2 * split] = FvBoundary::wall();
    faces[2 * split + 1] = FvBoundary::wall();
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-10;
    config.pressure_tolerance = 1e-12;
    config.pressure_profiles_pa[2 * axis] = Some(vec![
        offsets[0] + heads[0], offsets[0] + heads[0], 0.0,
        offsets[1] + heads[1], offsets[1] + heads[1],
    ]);
    config.pressure_profiles_pa[2 * axis + 1] = Some(vec![
        offsets[0], offsets[0], 0.0, offsets[1], offsets[1],
    ]);
    (domain, config)
}

/// G1/G3: disconnected Poiseuille flows depend on their own pressure heads,
/// including heads small compared with an unrelated component's datum.
#[test]
fn disconnected_pressure_components_preserve_small_heads_on_all_axes() {
    let gate = CancelGate::new();
    let offsets = [1_048_576.0, -1_099_511_627_776.0];
    for axis in 0..3 {
        let (domain, base_config) = separated_pressure_channels(axis, [0.0; 2], [0.5, -0.25]);
        let (_, shifted_config) = separated_pressure_channels(axis, offsets, [0.5, -0.25]);
        let base = simple_flow(&domain, &fluid(), &base_config, &gate).unwrap();
        let shifted = simple_flow(&domain, &fluid(), &shifted_config, &gate).unwrap();
        assert_eq!(shifted.velocity_m_s, base.velocity_m_s);
        assert_eq!(shifted.field, base.field);
        assert_eq!(shifted.report, base.report);
        let split = (axis + 1) % 3;
        let mut flow = [0.0; 2];
        for cell in 0..domain.cell_count() {
            if !domain.is_fluid(cell) {
                assert_eq!(shifted.pressure_pa[cell], 0.0);
                continue;
            }
            let at = domain.coords(cell);
            let region = usize::from(at[split] > 2);
            let pressure_error =
                shifted.pressure_pa[cell] - offsets[region] - base.pressure_pa[cell];
            assert!(pressure_error.abs() <= 2.0 * f64::EPSILON * offsets[region].abs());
            if domain.neighbor(cell, 2 * axis + 1).is_none() {
                flow[region] += shifted.field.outward(&domain, cell, 2 * axis + 1);
            }
        }
        // Each H=1, L=2, mu=1/4 duct has two transverse cells. The native
        // midpoint discretization gives Q=dp H^3/(12 mu L)*(1+2/2^2).
        assert!((flow[0] - 0.125).abs() < 1e-8, "{axis}: {flow:?}");
        assert!((flow[1] + 0.0625).abs() < 1e-8, "{axis}: {flow:?}");
    }
}

/// G0/G3: a pressure difference between disconnected regions does no work
/// and must not enter admission, momentum or the pressure residual.
#[test]
fn unequal_disconnected_constant_pressures_start_and_remain_at_rest() {
    let gate = CancelGate::new();
    for offsets in [[-1_048_576.0, 1_048_576.0], [-f64::MAX, f64::MAX]] {
        let (domain, config) = separated_pressure_channels(0, offsets, [0.0; 2]);
        let rest = simple_flow(&domain, &fluid(), &config, &gate).unwrap();
        assert!(rest.velocity_m_s.iter().all(|v| *v == [0.0; 3]));
        assert_eq!(rest.report.inflow_m3_s, 0.0);
        assert_eq!(rest.report.outflow_m3_s, 0.0);
        for cell in 0..domain.cell_count() {
            let at = domain.coords(cell);
            let expected = if at[1] == 2 { 0.0 } else { offsets[usize::from(at[1] > 2)] };
            assert_eq!(rest.pressure_pa[cell], expected);
        }
    }
}

/// G3: each source's pressure-flow work retains the same numerical datum as
/// its native momentum equations; large offsets cannot become viscous heat.
#[test]
fn disconnected_pressure_work_and_heat_ignore_independent_offsets() {
    use super::viscous_heat::{ViscousHeatControls, creeping_flow_with_pressure_work_heat};

    let (domain, mut config) = separated_pressure_channels(0, [0.0; 2], [0.5, -0.25]);
    config.tolerance = 1e-11;
    config.momentum_tolerance = 1e-10;
    let mu = vec![100.0; domain.cell_count()];
    let controls = ViscousHeatControls { maximum_reynolds: 0.1, power_tolerance_w: 1e-8 };
    let gate = CancelGate::new();
    let (flow, heat) =
        creeping_flow_with_pressure_work_heat(&domain, &fluid(), &mu, &config, controls, &gate)
            .unwrap();
    let mut supply = 0.0;
    for y in [0, 1] {
        supply += 0.5 * flow.field.flux_x(0, y, 0);
    }
    for y in [3, 4] {
        supply -= 0.25 * flow.field.flux_x(0, y, 0);
    }
    assert!(supply > 0.0);
    assert!((supply - heat.supplied_pressure_power_w).abs() < 1e-10);
    assert!((supply - heat.heat_w).abs() < controls.power_tolerance_w);
    for profile in config.pressure_profiles_pa.iter_mut().flatten() {
        for (row, pressure) in profile.iter_mut().enumerate() {
            if row < 2 {
                *pressure += 1_099_511_627_776.0;
            } else if row > 2 {
                *pressure -= 1_099_511_627_776.0;
            }
        }
    }
    let (shifted, shifted_heat) =
        creeping_flow_with_pressure_work_heat(&domain, &fluid(), &mu, &config, controls, &gate)
            .unwrap();
    assert_eq!(shifted.velocity_m_s, flow.velocity_m_s);
    assert_eq!(shifted_heat.cell_heat_w, heat.cell_heat_w);
    assert_eq!(shifted_heat.pressure_power_w, heat.pressure_power_w);
    assert!(shifted_heat.residual_w.abs() < controls.power_tolerance_w);
}

#[cfg(feature = "enthalpy-coupling")]
#[test]
fn component_pressure_checkpoint_retains_datums_and_rebases_joined_regions() {
    use crate::conjugate::Voxel;

    let (domain, config) = separated_pressure_channels(0, [1024.0, -2048.0], [0.5, -0.25]);
    let fluid = fluid();
    let gate = CancelGate::new();
    let mut solver = Solver::new(&domain, &fluid, &config);
    for _ in 0..5 {
        solver.sweep(&gate).unwrap();
    }
    let state = solver.momentum_state();
    let mut restored = Solver::new(&domain, &fluid, &config);
    restored.restore_momentum(&state).unwrap();
    assert_eq!(restored.pressure, solver.pressure);
    assert_eq!(restored.pressure_reference_pa, solver.pressure_reference_pa);
    assert_eq!(restored.vel, solver.vel);

    // Opening the original divider changes only previously fixed zero faces.
    // The new connected pressure field must use one common coordinate chart.
    let mut joined = domain.clone();
    for x in 0..8 {
        joined.set(x, 2, 0, Voxel::Fluid);
    }
    let mut restored = Solver::new(&joined, &fluid, &config);
    restored.restore_momentum(&state).unwrap();
    assert!(matches!(restored.pressure_reference_pa, PressureReferences::Uniform(_)));
    for cell in 0..domain.cell_count() {
        if domain.is_fluid(cell) {
            let before = state.pressure[cell] + state.pressure_reference_pa.at(cell);
            let after = restored.pressure[cell] + restored.pressure_reference_pa.at(cell);
            assert!((before - after).abs() <= 4.0 * f64::EPSILON * before.abs());
        }
    }
}

#[test]
fn fan_ambient_does_not_choose_an_unconnected_reservoir_datum() {
    use crate::conjugate::Voxel;

    let (mut domain, mut config) = separated_pressure_channels(0, [0.0, 1_048_576.0], [0.0; 2]);
    for y in [3, 4] {
        domain.set(0, y, 0, Voxel::Solid(0));
    }
    config.faces[0] = FvBoundary::Inlet { velocity: [0.1, 0.0, 0.0] };
    config.pressure_profiles_pa[0] = None;
    config.fan = Some(FanInlet {
        face: Face3::XMin,
        curve: FanCurve::new(&[(0.0, 1.0), (1.0, 0.0)]).unwrap(),
    });
    admit(&domain, &fluid(), &config).unwrap();
    let references = config.pressure_reference_pa(&domain);
    assert_eq!(references.at(domain.index(0, 0, 0)), 0.0);
    assert_eq!(references.at(domain.index(1, 3, 0)), 1_048_576.0);
}
