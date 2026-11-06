use super::*;
use crate::d3q27::D3q27;
use crate::nsf3::{NsfCharacteristicAxis, NsfCharacteristicBoundary, NsfSpec, WallThermal};
use fs_material::gas::GasSpec;

fn spec(dims: [usize; 3], periodic: [bool; 3]) -> NsfSpec {
    NsfSpec {
        dims,
        periodic,
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 300.0,
        pressure_pa: 100_000.0,
        dx_m: 0.001,
        temperature_window_k: [250.0, 400.0],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 0.3,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "Cramer 2012 N2 ratio; body-acceleration regression",
        maximum_cells: 2048,
        wall: WallThermal::Adiabatic,
        workers: 1,
        transport_points: 17,
    }
}

fn grid(dims: [usize; 3]) -> NsfGrid3<D3q27> {
    NsfGrid3::new(
        spec(dims, [true; 3]),
        |_| false,
        |_| (1.0, [0.01, -0.02, 0.015], 300.0),
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn uniform_body_drive_matches_gravity_and_does_mechanical_work_without_heat() {
    let acceleration = [9_000.0, -17_000.0, 12_000.0];
    let mut body = grid([3; 3]);
    let mut gravity = grid([3; 3]).with_gravity_si(acceleration).unwrap();
    let initial = body.cell([1; 3]);
    let mass = body.total_mass_kg();
    let energy = body.total_energy_j();
    let before = body.state_hash();
    body.set_body_acceleration_si(&vec![[0.0; 3]; body.n], || false)
        .unwrap();
    assert_eq!(before, body.state_hash());
    assert!(body.body_acceleration.is_empty());
    body.set_body_acceleration_si(&vec![acceleration; body.n], || false)
        .unwrap();
    assert_eq!(body.body_acceleration_m_s2([1; 3]).unwrap(), acceleration);
    assert_eq!(body.f, gravity.f);
    assert_eq!(body.force, gravity.force);
    near(body.total_energy_j(), energy, 0.0);
    let mut work = 0.0;
    for step in 1..=4 {
        let b = body.step().unwrap();
        let g = gravity.step().unwrap();
        assert_eq!(body.f, gravity.f);
        assert_eq!(body.e, gravity.e);
        assert_eq!(body.theta, gravity.theta);
        assert_eq!(b.body_force_work_j, g.gravity_work_j);
        let mut comparable = b;
        comparable.gravity_work_j = b.body_force_work_j;
        comparable.body_force_work_j = 0.0;
        assert_eq!(comparable, g);
        let (_, velocity, temperature) = body.cell([1; 3]);
        for axis in 0..3 {
            let change_m_s = (velocity[axis] - initial.1[axis]) * body.law.dx_m / body.law.dt_s;
            near(
                change_m_s,
                acceleration[axis] * (step as f64 * body.law.dt_s),
                5e-10,
            );
        }
        near(temperature, initial.2, 2e-8);
        near(body.total_mass_kg(), mass, 2e-13 * mass);
        work += b.body_force_work_j;
        near(body.total_energy_j() - energy, work, 2e-12 * energy);
        assert!(b.residual_rel < 2e-12);
        assert_eq!(b.source_energy_j, 0.0);
    }
    assert_eq!(body.cumulative_ledger().body_force_work_j, work);
}

#[test]
fn independent_gravity_body_and_manufactured_supplies_compose_once() {
    let initial = [0.01, -0.02, 0.015];
    let body_si = [8_000.0, -4_000.0, 3_000.0];
    let gravity_si = [-1_500.0, 5_000.0, -2_000.0];
    let manufactured = [2e-5, -1e-5, 1e-5];
    let reference = grid([3; 3]);
    let scale = reference.law.dt_s * (reference.law.dt_s / reference.law.dx_m);
    let body = body_si.map(|a| a * scale);
    let gravity = gravity_si.map(|a| a * scale);
    let total: [f64; 3] = core::array::from_fn(|a| body[a] + gravity[a] + manufactured[a]);
    let calls = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let observed = calls.clone();
    let source = Box::new(move |_: [usize; 3], time: f64| {
        observed.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        let power = (0..3)
            .map(|a| manufactured[a] * (initial[a] + total[a] * time))
            .sum();
        (0.0, manufactured, power)
    });
    let mut flow = NsfGrid3::<D3q27>::new_with_sources(
        spec([3; 3], [true; 3]),
        |_| false,
        |_| (1.0, initial, 300.0),
        Some(source),
    )
    .unwrap()
    .with_gravity_si(gravity_si)
    .unwrap();
    let accepted_calls = calls.load(std::sync::atomic::Ordering::SeqCst);
    flow.set_body_acceleration_si(&vec![body_si; flow.n], || false)
        .unwrap();
    assert_eq!(
        calls.load(std::sync::atomic::Ordering::SeqCst),
        accepted_calls
    );
    for step in 1..=4 {
        let ledger = flow.step().unwrap();
        let (_, velocity, temperature) = flow.cell([1; 3]);
        for axis in 0..3 {
            near(
                velocity[axis],
                initial[axis] + step as f64 * total[axis],
                3e-13,
            );
        }
        let expected_work = |acceleration: [f64; 3]| {
            flow.n as f64
                * flow.law.energy_unit_j()
                * (0..3)
                    .map(|a| acceleration[a] * (initial[a] + (step as f64 - 0.5) * total[a]))
                    .sum::<f64>()
        };
        for (actual, expected) in [
            (ledger.body_force_work_j, expected_work(body)),
            (ledger.gravity_work_j, expected_work(gravity)),
            (ledger.source_energy_j, expected_work(manufactured)),
        ] {
            near(actual, expected, 3e-11 * expected.abs().max(1e-20));
        }
        near(temperature, 300.0, 2e-8);
        assert!(ledger.residual_rel < 3e-12);
    }

    // Opposed mechanical controls do no NET work and must not become heat.
    let mut balanced = grid([3; 3]).with_gravity_si(gravity_si).unwrap();
    balanced
        .set_body_acceleration_si(&vec![gravity_si.map(|a| -a); balanced.n], || false)
        .unwrap();
    let energy = balanced.total_energy_j();
    for _ in 0..3 {
        let ledger = balanced.step().unwrap();
        near(ledger.gravity_work_j + ledger.body_force_work_j, 0.0, 1e-25);
        near(balanced.total_energy_j(), energy, 2e-12 * energy);
        near(balanced.cell([1; 3]).2, 300.0, 2e-8);
    }
}

#[test]
fn force_uses_evolving_compressible_density_with_mass_injection() {
    let mass_source = 0.001;
    let reference_internal = crate::nsf3::THETA0 / (grid([3; 3]).law.gamma - 1.0);
    let initial_velocity = [0.01, -0.02, 0.015];
    let mut flow = NsfGrid3::<D3q27>::new_with_sources(
        spec([3; 3], [true; 3]),
        |_| false,
        |_| (1.0, initial_velocity, 300.0),
        // Inject reference internal energy and no momentum. The native source
        // owner retains the resulting kinetic/internal exchange.
        Some(Box::new(move |_, _| {
            (mass_source, [0.0; 3], mass_source * reference_internal)
        })),
    )
    .unwrap();
    let acceleration_si = [8_000.0, -4_000.0, 3_000.0];
    let acceleration = acceleration_si.map(|a| a * flow.law.dt_s * (flow.law.dt_s / flow.law.dx_m));
    flow.set_body_acceleration_si(&vec![acceleration_si; flow.n], || false)
        .unwrap();
    for step in 1..=4 {
        let before = flow.cell([1; 3]);
        let ledger = flow.step().unwrap();
        let after = flow.cell([1; 3]);
        let time = step as f64;
        near(after.0, 1.0 + mass_source * time, 3e-14);
        for axis in 0..3 {
            let momentum = initial_velocity[axis]
                + acceleration[axis] * (time + 0.5 * mass_source * time * time);
            near(after.0 * after.1[axis], momentum, 3e-13);
        }
        let expected_work = 0.5
            * flow.n as f64
            * flow.law.energy_unit_j()
            * (0..3)
                .map(|a| acceleration[a] * (before.0 * before.1[a] + after.0 * after.1[a]))
                .sum::<f64>();
        near(
            ledger.body_force_work_j,
            expected_work,
            3e-11 * expected_work.abs(),
        );
        assert!(ledger.residual_rel < 3e-12);
    }
}

#[test]
fn a_spatial_drive_creates_the_resolved_vorticity_with_second_order_refinement() {
    let length = 0.012;
    let amplitude = 10_000.0;
    let mut errors = Vec::new();
    for n in [12, 24, 48] {
        let mut specification = spec([n, 3, 3], [true; 3]);
        specification.dx_m = length / n as f64;
        let mut flow =
            NsfGrid3::<D3q27>::new(specification, |_| false, |_| (1.0, [0.0; 3], 300.0)).unwrap();
        let k = core::f64::consts::TAU / n as f64;
        let drive: Vec<_> = (0..flow.n)
            .map(|i| {
                [
                    0.0,
                    amplitude * (k * (flow.at(i)[0] as f64 + 0.5)).sin(),
                    0.0,
                ]
            })
            .collect();
        flow.set_body_acceleration_si(&drive, || false).unwrap();
        for x in 0..n {
            for omega in flow.vorticity_budget([x, 1, 1]).unwrap().flow.vorticity() {
                near(omega, 0.0, 3e-16);
            }
        }
        let mass = flow.total_mass_kg();
        let ledger = flow.step().unwrap();
        let mut error = 0.0;
        let mut norm = 0.0;
        let source_scale = amplitude * core::f64::consts::TAU / length;
        for x in 0..n {
            let angle = k * (x as f64 + 0.5);
            let exact = source_scale * angle.cos();
            let centered = amplitude * k.sin() / flow.law.dx_m * angle.cos();
            let force_curl = flow.body_acceleration_curl_si([x, 1, 1]).unwrap();
            near(force_curl[2], centered, 2e-12 * source_scale);
            let q = flow.vorticity_budget([x, 1, 1]).unwrap();
            let rate = (q.flow.vorticity()[2] / flow.law.dt_s) / flow.law.dt_s;
            // Independent first-step D3Q27 transverse source filter:
            // half of a plus half of (2/3 + cos(k)/3)*a.
            let discrete = centered * (5.0 / 6.0 + k.cos() / 6.0);
            near(rate, discrete, 3e-8 * source_scale);
            error += (rate - exact).powi(2);
            norm += exact * exact;
        }
        errors.push((error / norm).sqrt());
        near(flow.total_mass_kg(), mass, 2e-13 * mass);
        assert!(ledger.body_force_work_j > 0.0);
        assert!(ledger.residual_rel < 3e-12);
    }
    assert!(errors[1] < 0.3 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.3 * errors[1], "{errors:?}");
    assert!(errors[2] < 0.006, "{errors:?}");
}

#[test]
fn rebind_retains_velocity_stress_energy_and_source_controls() {
    let mut flow = NsfGrid3::<D3q27>::new(
        spec([8, 3, 3], [true; 3]),
        |_| false,
        |x| {
            (
                1.0,
                [
                    0.0,
                    0.01 * (core::f64::consts::TAU * x[0] as f64 / 8.0).sin(),
                    0.005,
                ],
                300.0,
            )
        },
    )
    .unwrap()
    .with_gravity_si([0.0, -9.80665, 0.0])
    .unwrap();
    for _ in 0..3 {
        flow.step().unwrap();
    }
    let before: Vec<_> = (0..flow.n).map(|i| flow.cell(flow.at(i))).collect();
    let (rho, _) = flow.macro_fields(&flow.f);
    let stress = flow.stress_field(
        &flow.f,
        &flow.force,
        &rho,
        &flow.theta,
        &flow.material,
        &flow.sgs_transport,
    );
    assert!(stress.iter().flatten().any(|v| v.abs() > 1e-12));
    let energy = flow.e.clone();
    let temperature = flow.theta.clone();
    let ledger = flow.cumulative_ledger();
    let steps = flow.steps();
    let field: Vec<_> = (0..flow.n)
        .map(|i| [10.0 * i as f64, -20.0, 30.0])
        .collect();
    for field in [field.as_slice(), &[]] {
        flow.set_body_acceleration_si(field, || false).unwrap();
        let (rho, _) = flow.macro_fields(&flow.f);
        let after_stress = flow.stress_field(
            &flow.f,
            &flow.force,
            &rho,
            &flow.theta,
            &flow.material,
            &flow.sgs_transport,
        );
        for (before, after) in stress.iter().flatten().zip(after_stress.iter().flatten()) {
            near(*after, *before, 5e-14);
        }
        for (i, old) in before.iter().enumerate() {
            let new = flow.cell(flow.at(i));
            near(new.0, old.0, 5e-14);
            for a in 0..3 {
                near(new.1[a], old.1[a], 5e-14);
            }
        }
        assert_eq!(flow.e, energy);
        assert_eq!(flow.theta, temperature);
        assert_eq!(flow.cumulative_ledger(), ledger);
        assert_eq!(flow.steps(), steps);
        assert_eq!(flow.gravity_m_s2(), [0.0, -9.80665, 0.0]);
    }
}

#[test]
fn refusal_and_cancellation_preserve_the_entire_accepted_owner() {
    let dims = [9, 4, 4];
    let field = vec![[500.0, -700.0, 300.0]; dims.into_iter().product()];
    let mut counted = grid(dims);
    let mut polls = 0;
    counted
        .set_body_acceleration_si(&field, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls >= 6);
    for cancel_at in [1, polls / 2, polls] {
        let mut flow = grid(dims);
        let state = flow.state_hash();
        let force = flow.force.clone();
        let div = flow.div.clone();
        let ledger = flow.cumulative_ledger();
        let mut calls = 0;
        let result = flow.set_body_acceleration_si(&field, || {
            calls += 1;
            calls == cancel_at
        });
        assert!(result.is_err());
        assert_eq!(flow.state_hash(), state);
        assert_eq!(flow.force, force);
        assert_eq!(flow.div, div);
        assert_eq!(flow.cumulative_ledger(), ledger);
        assert!(flow.body_acceleration.is_empty());
        flow.set_body_acceleration_si(&field, || false).unwrap();
        assert_eq!(flow.state_hash(), counted.state_hash());
    }
    let mut flow = grid(dims);
    flow.set_body_acceleration_si(&field, || false).unwrap();
    let state = flow.state_hash();
    for invalid in [
        vec![[0.0; 3]; flow.n - 1],
        vec![[f64::NAN, 0.0, 0.0]; flow.n],
        vec![[f64::from_bits(1), 0.0, 0.0]; flow.n],
        vec![[f64::MAX, 0.0, 0.0]; flow.n],
    ] {
        assert!(flow.set_body_acceleration_si(&invalid, || false).is_err());
        assert_eq!(flow.state_hash(), state);
        assert_eq!(flow.body_acceleration_m_s2([1; 3]).unwrap(), field[0]);
    }
}

#[test]
fn acceleration_curl_excludes_density_gradients_and_refuses_missing_stencils() {
    let mut flow = NsfGrid3::<D3q27>::new(
        spec([5; 3], [false; 3]),
        |_| false,
        |x| (1.0 + 0.01 * x[0] as f64, [0.0; 3], 300.0),
    )
    .unwrap();
    flow.set_body_acceleration_si(&vec![[0.0, 20.0, 0.0]; flow.n], || false)
        .unwrap();
    assert_eq!(flow.body_acceleration_curl_si([2; 3]).unwrap(), [0.0; 3]);
    assert!(flow.body_acceleration_curl_si([0, 2, 2]).is_err());
    assert!(flow.body_acceleration_curl_si([5, 2, 2]).is_err());
    let mut solid = NsfGrid3::<D3q27>::new(
        spec([5; 3], [true; 3]),
        |x| x == [3, 2, 2],
        |_| (1.0, [0.0; 3], 300.0),
    )
    .unwrap();
    let mut field = vec![[0.0; 3]; solid.n];
    field[solid.index([3, 2, 2])] = [f64::NAN; 3];
    solid.set_body_acceleration_si(&field, || false).unwrap();
    assert!(solid.body_acceleration.is_empty());
    assert!(solid.body_acceleration_curl_si([2; 3]).is_err());
    assert!(solid.body_acceleration_m_s2([3, 2, 2]).is_err());
    let thin = grid([5, 2, 3]);
    assert!(thin.body_acceleration_curl_si([2, 1, 1]).is_err());
}

#[test]
fn characteristic_trials_retain_the_held_field_and_its_work() {
    let make = || {
        let mut flow = grid([12, 3, 3]);
        // Rebuild with a nonperiodic characteristic axis and matched initial
        // reference velocity; the transverse body drive varies along x.
        let velocity_unit = flow.law.dx_m / flow.law.dt_s;
        let face = NsfCharacteristicBoundary {
            reference_pressure_pa: 100_000.0,
            reference_temperature_k: 300.0,
            reference_velocity_m_s: [2.0, 0.0, 0.0],
            incoming_characteristic_m_s: 0.0,
        };
        flow = NsfGrid3::new(
            spec([12, 3, 3], [false, true, true]),
            |_| false,
            |_| (1.0, [2.0 / velocity_unit, 0.0, 0.0], 300.0),
        )
        .unwrap()
        .with_characteristic_axis(NsfCharacteristicAxis {
            axis: 0,
            faces: [face; 2],
            maximum_relative_perturbation: 0.05,
        })
        .unwrap();
        let field: Vec<_> = (0..flow.n)
            .map(|i| {
                [
                    100.0,
                    1000.0 * (core::f64::consts::TAU * (flow.at(i)[0] as f64 + 0.5) / 12.0).sin(),
                    0.0,
                ]
            })
            .collect();
        flow.set_body_acceleration_si(&field, || false).unwrap();
        flow
    };
    let mut direct = make();
    let mut previewed = make();
    let before = previewed.state_hash();
    direct.set_incoming_characteristics([-0.05, 0.02]).unwrap();
    let mut work = 0.0;
    for _ in 0..3 {
        work += direct.step().unwrap().body_force_work_j;
    }
    assert!(work > 0.0);
    {
        let trial = previewed
            .preview_characteristic_interval([-0.05, 0.02], 3, 324, || false)
            .unwrap();
        assert_eq!(trial.grid().state_hash(), direct.state_hash());
        assert_eq!(trial.energy_ledger().body_force_work_j, work);
        assert!(trial.energy_ledger().residual_rel < 3e-12);
    }
    assert_eq!(previewed.state_hash(), before);
    let trial = previewed
        .preview_characteristic_interval([-0.05, 0.02], 3, 324, || false)
        .unwrap();
    trial.accept();
    assert_eq!(previewed.state_hash(), direct.state_hash());
    assert_eq!(previewed.cumulative_ledger(), direct.cumulative_ledger());
}

mod unsteady;
