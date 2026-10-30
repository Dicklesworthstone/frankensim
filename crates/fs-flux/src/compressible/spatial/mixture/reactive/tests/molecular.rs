use super::*;
use crate::compressible::spatial::mixture::molecular::{
    DissipationControl, MolecularTransport3d, SpatialTransport,
};
use fs_gas::multicomponent::dissipation::{
    ConstantTransport, TransportCoefficients, TransportDomain,
};
use fs_qty::DynViscosity;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;

fn operator(flow: &MixtureEuler3d, diffusion: f64) -> MolecularTransport3d<ConstantTransport> {
    let law = ConstantTransport::new(
        flow.gas().source().clone(),
        "manufactured-spatial-reacting-transport",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 0.001,
            maximum_pressure_pa: 100.0,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(0.03),
            bulk_viscosity: DynViscosity::new(0.01),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.02),
        },
        vec![0.0, diffusion, diffusion, 0.0],
    )
    .unwrap();
    let source = SpatialTransport::new(flow.gas().clone(), law, || false).unwrap();
    let bytes = MolecularTransport3d::required_workspace_bytes(flow, &source).unwrap();
    MolecularTransport3d::new(flow, source, bytes, || false).unwrap()
}

fn budget(flow: &ReactingEuler3d, molecular: &MolecularTransport3d<ConstantTransport>) -> usize {
    flow.workspace_bytes() + molecular.workspace_bytes()
}

fn loads(flow: &ReactingEuler3d) -> Vec<[f64; 6]> {
    BoundarySide::ALL
        .into_iter()
        .flat_map(|side| {
            flow.boundary_impulses(side).map(|q| {
                [
                    q.species_masses[0],
                    q.species_masses[1],
                    q.momentum[0],
                    q.momentum[1],
                    q.momentum[2],
                    q.energy,
                ]
            })
        })
        .collect()
}

fn molecular_loads(molecular: &MolecularTransport3d<ConstantTransport>) -> Vec<[f64; 6]> {
    BoundarySide::ALL
        .into_iter()
        .flat_map(|side| {
            molecular.boundary_impulses(side).map(|q| {
                [
                    q.species_masses[0],
                    q.species_masses[1],
                    q.momentum[0],
                    q.momentum[1],
                    q.momentum[2],
                    q.energy,
                ]
            })
        })
        .collect()
}

#[test]
fn g3_uniform_reacting_molecular_field_retains_chemistry_and_heat_driven_cfl_restarts() {
    let mut reference = uniform();
    let expected = reference
        .step(Time::new(0.2), ReactiveStepControl::default(), || false)
        .unwrap();
    let mut flow = uniform();
    let mut molecular = operator(flow.flow(), 0.02);
    let bytes = budget(&flow, &molecular);
    let report = flow
        .step_with_dissipation(
            &mut molecular,
            Time::new(0.2),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            || false,
        )
        .unwrap();
    balance(&report);
    assert!(report.coupling_retries > 0);
    assert_eq!(report.elapsed, expected.elapsed);
    assert_eq!(report.transport, expected.transport);
    assert_eq!(report.chemistry, expected.chemistry);
    assert_eq!(flow.flow().states(), reference.flow().states());
    let split = report.dissipation.unwrap();
    assert_eq!(split.first.elapsed.value(), 0.5 * report.elapsed.value());
    assert_eq!(split.second.elapsed, split.first.elapsed);
    assert_eq!(molecular.last_step(), Some(&split.second));
    assert_eq!(budget(&flow, &molecular), bytes);
}

fn reacting_diffusion(n: usize) -> f64 {
    let (gas, source) = mechanism(0.0, 2.0, 0.0);
    let dx = 1.0 / n as f64;
    let z = std::f64::consts::PI * dx;
    let average = (z.sin() / z).powi(3);
    let diffusion = 0.02;
    let velocity = [0.4, 0.25, -0.15];
    // Original reversible A <-> B source gives d(Y_A-1/2)/dt=-4(Y_A-1/2).
    // Three equal wave numbers add 12*pi^2*D to that exact decay rate.
    let decay_rate = 4.0 + diffusion * 12.0 * std::f64::consts::PI.powi(2);
    let exact = |x: [f64; 3], time: f64| {
        0.5 + 0.2
            * average
            * (-decay_rate * time).exp()
            * (2.0 * std::f64::consts::PI * (x[0] + x[1] + x[2] - 0.5 * time)).sin()
    };
    let grid = CartesianGrid::new(n, n, n, dx, dx, dx).unwrap();
    let mut flow = coupled(
        field(gas, grid, PERIODIC, |x| {
            let a = exact(x, 0.0);
            ([a, 1.0 - a], velocity, 1.0)
        }),
        source,
    );
    let mut molecular = operator(flow.flow(), diffusion);
    let bytes = budget(&flow, &molecular);
    let end = 0.03;
    let mut steps = 0;
    while flow.flow().time().value() < end {
        let report = flow
            .step_with_dissipation(
                &mut molecular,
                Time::new(end - flow.flow().time().value()),
                ReactiveStepControl::default(),
                DissipationControl::default(),
                bytes,
                || false,
            )
            .unwrap();
        balance(&report);
        let split = report.dissipation.unwrap();
        assert!(split.first.maximum_parabolic_number <= 0.40000000000001);
        assert!(split.second.maximum_parabolic_number <= 0.40000000000001);
        steps += 1;
        assert!(steps < 100);
    }
    let mut error = 0.0;
    for i in 0..flow.flow().cell_count() {
        let q = flow.flow().cell(i).unwrap();
        let p = q.primitive(flow.flow().gas(), || false).unwrap();
        close(p.pressure, 1.0, 2e-10);
        close(p.temperature, 1.0, 2e-10);
        for axis in 0..3 {
            close(p.velocity[axis], velocity[axis], 2e-10);
        }
        error += (q.partial_densities[0] - exact(grid.cell_center(i).unwrap(), end)).abs();
    }
    error / flow.flow().cell_count() as f64
}

#[test]
fn g1_oblique_three_dimensional_reaction_advection_diffusion_refines_to_exact_solution() {
    let coarse = reacting_diffusion(12);
    let fine = reacting_diffusion(24);
    let order = (coarse / fine).ln() / 2.0_f64.ln();
    eprintln!("3-D reacting diffusion L1: n12={coarse:e}, n24={fine:e}, order={order:e}");
    assert!(coarse < 0.01);
    assert!(
        order > 1.8,
        "coarse={coarse:e}, fine={fine:e}, order={order:e}"
    );
}

#[test]
fn g4_all_molecular_reacting_phases_and_support_loads_publish_atomically() {
    let create = || {
        let (gas, source) = mechanism(5.0, 8.0, 0.0);
        let grid = CartesianGrid::new(3, 3, 3, 0.2, 0.25, 0.3).unwrap();
        let input = field(
            gas,
            grid,
            Boundaries3d {
                x: WALLS,
                y: WALLS,
                z: WALLS,
            },
            |x| {
                let a = 0.5 + 0.2 * x[0] - 0.1 * x[2];
                (
                    [a, 1.0 - a],
                    [0.2 + 0.1 * x[1], -0.15, 0.12],
                    1.0 + 0.1 * x[2],
                )
            },
        );
        // An oversized incoming buffer must remain counted as it moves between
        // the two flow fields and the molecular operator during accepted swaps.
        let mut states = Vec::with_capacity(input.states.len() + 31);
        states.extend_from_slice(&input.states);
        let input = MixtureEuler3d::new(
            input.gas().clone(),
            *input.grid(),
            states,
            input.boundaries(),
            crate::compressible::spatial::mixture::MAX_MIXTURE_WORKSPACE_BYTES,
            || false,
        )
        .unwrap();
        let flow = coupled(input, source);
        let molecular = operator(flow.flow(), 0.02);
        (flow, molecular)
    };
    let (mut expected, mut expected_molecular) = create();
    let (mut flow, mut molecular) = create();
    let bytes = budget(&flow, &molecular);
    let controls = DissipationControl {
        maximum_step_s: 0.0005,
        ..DissipationControl::default()
    };
    for (flow, molecular) in [
        (&mut expected, &mut expected_molecular),
        (&mut flow, &mut molecular),
    ] {
        flow.step_with_dissipation(
            molecular,
            Time::new(0.001),
            ReactiveStepControl::default(),
            controls,
            bytes,
            || false,
        )
        .unwrap();
        assert_eq!(budget(flow, molecular), bytes);
    }
    let mut polls = 0;
    let reference = expected
        .step_with_dissipation(
            &mut expected_molecular,
            Time::new(0.002),
            ReactiveStepControl::default(),
            controls,
            bytes,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    balance(&reference);
    let split = reference.dissipation.unwrap();
    assert_ne!(
        reference.boundary_change.momentum,
        reference.transport.boundary_change.momentum
    );
    let mut outward = [0.0; 3];
    for load in loads(&expected) {
        assert_eq!(load[0], 0.0);
        assert_eq!(load[1], 0.0);
        assert_eq!(load[5], 0.0);
        for axis in 0..3 {
            outward[axis] += load[axis + 2];
        }
    }
    for axis in 0..3 {
        close(
            reference.after.momentum[axis] - reference.before.momentum[axis] + outward[axis],
            0.0,
            3e-12,
        );
        close(
            reference.boundary_change.momentum[axis],
            reference.transport.boundary_change.momentum[axis]
                + split.first.boundary_change.momentum[axis]
                + split.second.boundary_change.momentum[axis],
            3e-12,
        );
    }
    let states = flow.flow().states().to_vec();
    let time = flow.flow().time();
    let report = *flow.last_step().unwrap();
    let wave_report = *flow.flow().last_step().unwrap();
    let molecular_report = *molecular.last_step().unwrap();
    let accepted_loads = loads(&flow);
    let accepted_molecular_loads = molecular_loads(&molecular);
    let mut count = 0;
    assert_eq!(
        flow.step_with_dissipation(
            &mut molecular,
            Time::new(0.002),
            ReactiveStepControl::default(),
            controls,
            bytes,
            || {
                count += 1;
                count >= polls
            },
        ),
        Err(ReactiveError::Cancelled)
    );
    for refused in [
        flow.step_with_dissipation(
            &mut molecular,
            Time::new(0.001),
            ReactiveStepControl {
                maximum_coupling_retries: 0,
                ..ReactiveStepControl::default()
            },
            DissipationControl {
                maximum_step_s: 1e-6,
                maximum_steps: 1,
                ..controls
            },
            bytes,
            || false,
        ),
        flow.step_with_dissipation(
            &mut molecular,
            Time::new(0.001),
            ReactiveStepControl::default(),
            controls,
            bytes - 1,
            || false,
        ),
    ] {
        assert!(refused.is_err());
    }
    assert_eq!(flow.flow().states(), states);
    assert_eq!(flow.flow().time(), time);
    assert_eq!(flow.last_step(), Some(&report));
    assert_eq!(flow.flow().last_step(), Some(&wave_report));
    assert_eq!(molecular.last_step(), Some(&molecular_report));
    assert_eq!(loads(&flow), accepted_loads);
    assert_eq!(molecular_loads(&molecular), accepted_molecular_loads);
    assert_eq!(budget(&flow, &molecular), bytes);
    let replay = flow
        .step_with_dissipation(
            &mut molecular,
            Time::new(0.002),
            ReactiveStepControl::default(),
            controls,
            bytes,
            || false,
        )
        .unwrap();
    assert_eq!(replay, reference);
    assert_eq!(flow.flow().states(), expected.flow().states());
    assert_eq!(loads(&flow), loads(&expected));
    assert_eq!(
        molecular_loads(&molecular),
        molecular_loads(&expected_molecular)
    );
    assert_eq!(budget(&flow, &molecular), bytes);
}
