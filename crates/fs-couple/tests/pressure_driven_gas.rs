//! Native joint mechanics/gas proof for nonlinear pressure supply work.
//! The original fs-tribo storage and transport run inside the existing fs-phs
//! implicit equation, with no separate gas clock or endpoint pressure force.
use fs_phs::{NonlinearPortPower, PortHamiltonian, PreparedStepError, StepWorkspace, Storage};
use fs_tribo::resistive_film::isothermal::reservoirs::{GasChannelFlow, PressureDrivenGasFilm};
use fs_tribo::resistive_film::isothermal::{AmbientGas, IsothermalFilm};
use fs_tribo::resistive_film::{FilmCell, FilmChannel, FilmLimits, GapPort, ResistiveFilm};
use std::rc::Rc;

const MASS: f64 = 0.2;
const RT: f64 = 293.15 * 287.05;
struct JointStorage(Rc<PressureDrivenGasFilm>);
impl Storage for JointStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.0.observe(&[x[0]], &x[2..]).map_or(f64::NAN, |gas| {
            gas.free_energy_j + x[1] * x[1] / (2.0 * MASS)
        })
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        let mut q = [0.0];
        let mut z = [0.0];
        if self
            .0
            .gradient_into(&[x[0]], &x[2..], &mut q, &mut z)
            .is_err()
        {
            out.fill(f64::NAN);
            return;
        }
        out.copy_from_slice(&[q[0], x[1] / MASS, z[0]]);
    }
}
fn model() -> (Rc<PressureDrivenGasFilm>, PortHamiltonian) {
    let gap = GapPort {
        reference_m: 5e-5,
        closure: vec![1.0],
    };
    let geometry = ResistiveFilm::new(
        vec![FilmCell {
            area_m2: 0.002,
            gap: gap.clone(),
        }],
        vec![FilmChannel {
            from: 0,
            to: None,
            width_m: 0.003,
            length_m: 0.01,
            gap,
        }],
        1,
        1.8e-5,
        FilmLimits {
            minimum_cell_gap_m: 1e-6,
            maximum_gap_m: 0.01,
            maximum_pressure_pa: 300_000.0,
        },
    )
    .unwrap();
    let ambient = IsothermalFilm::at_ambient(
        geometry,
        &[0.0],
        AmbientGas {
            pressure_pa: 100_000.0,
            temperature_k: 293.15,
            specific_gas_constant_j_kg_k: 287.05,
        },
    )
    .unwrap();
    let film = Rc::new(PressureDrivenGasFilm::new(ambient, &[Some(150_000.0)]).unwrap());
    let system = PortHamiltonian::new(
        3,
        0,
        vec![0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        vec![0.0; 9],
        vec![],
        Box::new(JointStorage(Rc::clone(&film))),
    )
    .unwrap();
    (film, system)
}

#[test]
fn supplied_gas_mass_and_pressure_motion_share_one_discrete_energy_step() {
    let (film, system) = model();
    let port = |x: &[f64], e: &[f64], out: &mut [f64]| {
        out.fill(0.0);
        film.flow_into(
            &[x[0]],
            &x[2..],
            &e[2..],
            &mut out[2..],
            &mut [GasChannelFlow::default()],
        )
        .ok()
        .map(|r| NonlinearPortPower {
            supplied_power_w: r.reservoir_power_w,
            dissipated_power_w: r.dissipated_power_w,
        })
    };
    let hessian = |x: &[f64], d: &[f64], out: &mut [f64]| {
        let (mut q, mut z) = ([0.0], [0.0]);
        if film
            .hessian_into(&[x[0]], &x[2..], &[d[0]], &d[2..], &mut q, &mut z)
            .is_err()
        {
            return false;
        }
        out.copy_from_slice(&[q[0], d[1] / MASS, z[0]]);
        true
    };
    let tangent = |x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]| {
        out.fill(0.0);
        film.flow_tangent_into(
            &[x[0]],
            &x[2..],
            &e[2..],
            &[dx[0]],
            &dx[2..],
            &de[2..],
            &[Some(0.0)],
            &mut out[2..],
        )
        .is_ok()
    };
    let mut analytic = StepWorkspace::new(&system).unwrap();
    let mut fd = StepWorkspace::new(&system).unwrap();
    let mut x = [0.0, 0.0, film.initial_state()[0]];
    let mut reference = x;
    let initial = film.observe(&[x[0]], &x[2..]).unwrap();
    let (mut supplied, mut dissipated, mut mass_in) = (0.0, 0.0, 0.0);
    let dt = 1e-5;
    for _ in 0..100 {
        let mut next = [0.0; 3];
        let mut expected = [0.0; 3];
        let receipt = analytic
            .step_into_supplied_controlled(
                &system,
                &x,
                &[],
                dt,
                &mut next,
                &mut [],
                &port,
                Some(&hessian),
                Some(&tangent),
                || false,
            )
            .unwrap();
        fd.step_into_supplied_controlled(
            &system,
            &reference,
            &[],
            dt,
            &mut expected,
            &mut [],
            &port,
            None,
            None,
            || false,
        )
        .unwrap();
        for i in 0..3 {
            assert!((next[i] - expected[i]).abs() < 2e-8);
        }
        let mid = core::array::from_fn::<_, 3, _>(|i| f64::midpoint(x[i], next[i]));
        let effort = analytic.last_step_discrete_effort().unwrap();
        let report = film
            .flow_into(
                &[mid[0]],
                &mid[2..],
                &effort[2..],
                &mut [0.0],
                &mut [GasChannelFlow::default()],
            )
            .unwrap();
        assert!((receipt.supplied - dt * report.reservoir_power_w).abs() < 2e-13);
        assert!((receipt.dissipated - dt * report.dissipated_power_w).abs() < 2e-13);
        assert!(receipt.dissipated >= 0.0);
        assert!(receipt.balance_residual().abs() < 2e-10);
        supplied += receipt.supplied;
        dissipated += receipt.dissipated;
        mass_in += dt * report.reservoir_mass_rate_kg_s;
        x = next;
        reference = expected;
    }
    let final_gas = film.observe(&[x[0]], &x[2..]).unwrap();
    assert!(
        x[0] < -1e-6 && x[1] < 0.0,
        "pressure must open the moving gap: {x:?}"
    );
    assert!(final_gas.mass_kg > initial.mass_kg);
    assert!(supplied > 0.0 && dissipated > 0.0);
    assert!((final_gas.mass_kg - initial.mass_kg - mass_in).abs() < 2e-13);
    assert!((system.hamiltonian(&x) - supplied + dissipated).abs() < 2e-8);
    // Independent conversion of the ORIGINAL scaled state to physical mass.
    assert!((final_gas.mass_kg - film.initial_state()[0] * x[2] / RT).abs() < 1e-20);
}

#[test]
fn pressure_supply_refusal_and_final_cancel_preserve_joint_outputs_for_retry() {
    let (film, system) = model();
    let port = |x: &[f64], e: &[f64], out: &mut [f64]| {
        out.fill(0.0);
        film.flow_into(
            &[x[0]],
            &x[2..],
            &e[2..],
            &mut out[2..],
            &mut [GasChannelFlow::default()],
        )
        .ok()
        .map(|r| NonlinearPortPower {
            supplied_power_w: r.reservoir_power_w,
            dissipated_power_w: r.dissipated_power_w,
        })
    };
    let x = [0.0, 0.0, film.initial_state()[0]];
    let mut clean = StepWorkspace::new(&system).unwrap();
    let mut expected = [0.0; 3];
    let mut polls = 0;
    let receipt = clean
        .step_into_supplied_controlled(
            &system,
            &x,
            &[],
            1e-5,
            &mut expected,
            &mut [],
            &port,
            None,
            None,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut work = StepWorkspace::new(&system).unwrap();
    let mut output = [19.0; 3];
    let mut seen = 0;
    assert_eq!(
        work.step_into_supplied_controlled(
            &system,
            &x,
            &[],
            1e-5,
            &mut output,
            &mut [],
            &port,
            None,
            None,
            || {
                seen += 1;
                seen == polls
            }
        ),
        Err(PreparedStepError::Cancelled)
    );
    assert_eq!(output, [19.0; 3]);
    let mut invalid = x;
    invalid[2] = 0.0;
    assert!(
        work.step_into_supplied_controlled(
            &system,
            &invalid,
            &[],
            1e-5,
            &mut output,
            &mut [],
            &port,
            None,
            None,
            || false
        )
        .is_err()
    );
    assert_eq!(output, [19.0; 3]);
    let retry = work
        .step_into_supplied_controlled(
            &system,
            &x,
            &[],
            1e-5,
            &mut output,
            &mut [],
            &port,
            None,
            None,
            || false,
        )
        .unwrap();
    assert_eq!(output, expected);
    assert_eq!(retry, receipt);
}
