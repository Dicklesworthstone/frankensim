//! Exact terminal trajectories retain the original rational coordinates.
use super::*;
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, PortId};
use fs_duct::{Duct, Segment};
use fs_material::gas::GasSpec;
use fs_vfit::mimo::{
    MimoModel,
    load::{MimoLoadLimits, PreparedMimoLoad},
};
use fs_vfit::model::PoleTerm;

#[derive(Clone, Copy)]
struct Resonator {
    r: f64,
    k: f64,
    a: f64,
    omega: f64,
    beta: f64,
}
fn source() -> (GasState, Resonator) {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * 0.001_f64.powi(2));
    let (omega, a) = (1300.0_f64, 25.0_f64);
    (
        gas,
        Resonator {
            r: 0.02 * zc,
            k: 0.2 * zc * omega,
            a,
            omega,
            beta: det::sqrt(omega * omega - a * a),
        },
    )
}
fn terminated_load(rate: u32, walls: bool, r: f64) -> DgPortLoad {
    let (gas, source) = source();
    // Independent R + K*s/(s^2 + 2*a*s + omega^2), not a fitted terminal.
    let model = MimoModel {
        ports: 1,
        poles: vec![PoleTerm::Pair {
            pole: C64::new(-source.a, source.beta),
            residue: C64::ZERO,
        }],
        residues: vec![vec![C64::new(
            0.5 * source.k,
            0.5 * source.k * source.a / source.beta,
        )]],
        d: vec![r],
        e: vec![0.0],
        reciprocal: true,
    };
    let dt = 1.0 / f64::from(rate);
    let terminal = PreparedMimoLoad::new(
        &model,
        dt,
        MimoLoadLimits {
            band_rad_s: [1.0, core::f64::consts::PI / dt],
            sampled_passivity_points: 129,
            max_states: 2,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap();
    let duct = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.001,
                length: 0.05,
            }],
        },
        &gas,
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Port(PortId(7)),
        DgBoundary::Port(PortId(9)),
    )
    .unwrap();
    let duct = if walls {
        duct.with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap()
    } else {
        duct
    };
    DgPortLoad::new_with_terminal(duct, dt, 1e-10, 256, terminal).unwrap()
}
fn terminated(rate: u32, memory: bool) -> SharedCavityDg {
    cavity(rate, memory, 1e6)
        .with_exact_downstream_duct(
            terminated_load(rate, memory, source().1.r),
            budget(),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
}
fn same_terminal(a: &SharedCavityDg, b: &SharedCavityDg) {
    same_state(a, b);
    for (x, y) in a
        .load()
        .terminal()
        .unwrap()
        .pole_states()
        .iter()
        .zip(b.load().terminal().unwrap().pole_states())
    {
        close(x.re, y.re, 1e-20);
        close(x.im, y.im, 1e-20);
    }
}

#[test]
fn full_terminal_trajectory_and_actual_port_work_survive_three_clock_partitions() {
    let gate = CancelGate::new_clock_free();
    let mut whole = terminated(24000, false);
    let mut halves = terminated(48000, false);
    let mut quarters = terminated(96000, false);
    let mut covariance = 0.0_f64;
    let mut terminal_work = 0.0_f64;
    for _ in 0..24 {
        let frame = whole.step(&[0.001], 0.1, &gate).unwrap().clone();
        let original = *whole.load().last_terminal_frame().unwrap();
        for (system, count) in [(&mut halves, 2), (&mut quarters, 4)] {
            let mut work = 0.0;
            for _ in 0..count {
                work += system.step(&[0.001], 0.1, &gate).unwrap().terminal_work_j;
            }
            same_terminal(&whole, system);
            close(frame.terminal_work_j, work, 2e-15);
        }
        close(original.port_work_j, original.load_work_j, 1e-20);
        close(frame.terminal_work_j, original.load_work_j, 1e-20);
        covariance = covariance.max(
            (original.port_work_j
                - original.outward_flow_m3_s * original.mean_pressure_pa / 24000.0)
                .abs(),
        );
        terminal_work += frame.terminal_work_j.abs();
        assert!(frame.balance_residual_j.abs() <= frame.balance_tolerance_j);
    }
    assert!(
        terminal_work > 0.0 && covariance > 1e-18,
        "{terminal_work:e} / {covariance:e}"
    );
    assert!(whole.load().terminal().unwrap().pole_states()[0].abs() > 0.0);
}

// The independent physical storage is specific to this authored RLC source;
// arbitrary fitted terminal coordinates are deliberately not assigned energy.
fn rlc_energy_and_heat(system: &SharedCavityDg) -> (f64, f64) {
    let (_, rlc) = source();
    let load = system.load();
    let z = load.terminal().unwrap().pole_states()[0];
    let q = z.im / rlc.beta;
    let v = z.re - rlc.a * q;
    let original = load.duct().continuous_operator();
    let (_, pressure) = original
        .ports()
        .iter()
        .find(|(id, _)| *id == PortId(9))
        .unwrap();
    let endpoint = original.energy_coordinates(load.state()).unwrap();
    let p: f64 = pressure.iter().zip(endpoint).map(|(b, x)| b * x).sum();
    let effort = rlc.k * v;
    let flow = (p - effort) / rlc.r;
    (
        0.5 * rlc.k * (v * v + rlc.omega.powi(2) * q * q),
        rlc.r * flow * flow + 2.0 * rlc.a * rlc.k * v * v,
    )
}

#[test]
fn original_rlc_storage_and_heat_independently_close_work_with_real_wall_memory() {
    let gate = CancelGate::new_clock_free();
    let rate = 24000;
    let subdivisions = 256;
    let mut coarse = terminated(rate, true);
    let mut fine = terminated(rate * subdivisions, true);
    let initial = coarse.stored_energy_j().unwrap();
    let limits = *coarse.load().terminal().unwrap().limits();
    let feedback_rate = {
        let original = coarse.load().duct().continuous_operator();
        let (_, column) = original
            .ports()
            .iter()
            .find(|(id, _)| *id == PortId(9))
            .unwrap();
        column.iter().map(|b| b * b).sum::<f64>() / source().1.r
    };
    eprintln!(
        "original outlet decay rate {feedback_rate:e}/s; lambda*h at 64/128/256: {:?}",
        [64, 128, 256].map(|panels| feedback_rate / f64::from(rate * panels))
    );
    let (mut authored, mut physical_loss, mut terminal_work) = (0.0, 0.0, 0.0);
    let mut heats = [0.0; 3];
    let mut fine_terminal_work = 0.0;
    for _ in 0..16 {
        let frame = coarse.step(&[0.001], 0.03, &gate).unwrap().clone();
        // The measured outlet rate is about 3.1e6/s. Simpson16/32/64 has
        // lambda*h about 8/4/2; its startup-heat error crosses zero and is not
        // asymptotic. Nested Simpson64/128/256 resolves that fast boundary layer
        // along ONE exact trajectory, without changing the physical tolerance.
        let mut quadrature = [rlc_energy_and_heat(&fine).1; 3];
        for step in 1..=subdivisions {
            fine_terminal_work += fine.step(&[0.001], 0.03, &gate).unwrap().terminal_work_j;
            let power = rlc_energy_and_heat(&fine).1;
            for (rung, stride) in [4, 2, 1].into_iter().enumerate() {
                if step % stride == 0 {
                    let index = step / stride;
                    let weight = if step == subdivisions {
                        1.0
                    } else if index % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    quadrature[rung] += weight * power;
                }
            }
        }
        for (rung, stride) in [4, 2, 1].into_iter().enumerate() {
            heats[rung] +=
                quadrature[rung] * f64::from(stride) / (3.0 * f64::from(rate * subdivisions));
        }
        same_terminal(&coarse, &fine);
        authored += frame.authored_external_work_j;
        physical_loss += frame.inferred_total_dissipation_j;
        terminal_work += frame.terminal_work_j;
        assert!(frame.modal.material_dissipation_j > 0.0);
    }
    let stored_terminal = rlc_energy_and_heat(&coarse).0;
    let errors = heats.map(|heat| (terminal_work - stored_terminal - heat).abs());
    let signed = heats.map(|heat| terminal_work - stored_terminal - heat);
    eprintln!(
        "terminal heat at 64/128/256: {heats:?}; signed residuals {signed:?}; coarse work {terminal_work:e}, fine work {fine_terminal_work:e}, difference {:e}; coarse/fine storage {stored_terminal:e}/{:e}",
        terminal_work - fine_terminal_work,
        rlc_energy_and_heat(&fine).0
    );
    assert!(errors[1] < errors[0] && errors[2] < errors[1], "{errors:?}");
    assert!((terminal_work - fine_terminal_work).abs() < 1e-18);
    let heat = heats[2];
    let tolerance = 1e-17 + 2e-6 * terminal_work.abs().max(stored_terminal + heat);
    assert!(
        (terminal_work - stored_terminal - heat).abs() < tolerance,
        "terminal work={terminal_work:e}, physical storage={stored_terminal:e}, heat={heat:e}"
    );
    assert!(
        (coarse.stored_energy_j().unwrap() + stored_terminal + physical_loss + heat
            - initial
            - authored)
            .abs()
            < 2e-14 + tolerance
    );
    assert!(stored_terminal > 0.0 && heat > 0.0);
    assert!(coarse.load().state().wall_stored_energy_j() > 0.0);
    assert_eq!(*coarse.load().terminal().unwrap().limits(), limits);
    // A failed late candidate preserves all three original accepted owners.
    let frame = coarse.last_frame().unwrap().clone();
    let duct = coarse.load().state().clone();
    let poles = coarse.load().terminal().unwrap().pole_states().to_vec();
    let state = endpoint(&coarse);
    let memory = coarse.cavity().network().material_memory().to_vec();
    assert!(coarse.step(&[1e20], 0.0, &gate).is_err());
    let cancel = CancelGate::new_clock_free();
    cancel.request();
    assert!(coarse.step(&[0.0], 0.0, &cancel).is_err());
    assert_eq!(coarse.last_frame().unwrap(), &frame);
    assert_eq!(coarse.load().state(), &duct);
    assert_eq!(coarse.load().terminal().unwrap().pole_states(), poles);
    assert_eq!(endpoint(&coarse), state);
    assert_eq!(coarse.cavity().network().material_memory(), memory);
}

#[test]
fn strictly_proper_terminal_refuses_instantaneous_elimination_without_invented_storage() {
    let gate = CancelGate::new_clock_free();
    // Original held-flow mean coupling has a valid nonzero dynamic impedance,
    // but its zero instantaneous D requires a descriptor-system treatment.
    let rate = 48000;
    let dt = 1.0 / f64::from(rate);
    let (gas, _) = source();
    let zc = gas.density * gas.sound_speed / (core::f64::consts::PI * 0.001_f64.powi(2));
    // Same admitted strictly proper primitive as the harmonic owner test:
    // Z(s)=Zc*3000/(s+3000), with a genuinely zero instantaneous D.
    let source = MimoModel {
        ports: 1,
        poles: vec![PoleTerm::Real {
            pole: -3000.0,
            residue: 0.0,
        }],
        residues: vec![vec![C64::from_re(zc * 3000.0)]],
        d: vec![0.0],
        e: vec![0.0],
        reciprocal: true,
    };
    let terminal = PreparedMimoLoad::new(
        &source,
        dt,
        MimoLoadLimits {
            band_rad_s: [100.0, core::f64::consts::PI / dt],
            sampled_passivity_points: 33,
            max_states: 8,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap();
    let original = terminated_load(rate, false, 0.02 * zc);
    let load =
        DgPortLoad::new_with_terminal(original.duct().clone(), dt, 1e-10, 256, terminal).unwrap();
    let mean = cavity(rate, false, 1e6)
        .with_downstream_duct(load, &gate)
        .unwrap();
    assert!(mean.with_exact_evolution(budget(), &gate).is_err());
}
