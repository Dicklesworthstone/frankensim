use super::*;
use crate::compressible::mixture::tests::{close, gas as make_gas, source};
use crate::compressible::spatial::mixture::state_into;
use crate::compressible::spatial::{
    AxisBoundaries, Boundaries3d, BoundaryCondition, CartesianGrid,
};
use fs_gas::chemistry::{IsochoricKineticsMethod, advance_adiabatic};
use fs_gas::multicomponent::MixtureGas;
use fs_qty::{Amount, Temperature};
use fs_thermochem::kinetics::{
    ArrheniusRateV1, ElementaryReactionV1, ReactionProgressRateQuantityV1,
};
use fs_thermochem::{
    MolarEnergyQuantityV1, ReactionId, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

const PERIODIC: Boundaries3d = Boundaries3d {
    x: AxisBoundaries::Periodic,
    y: AxisBoundaries::Periodic,
    z: AxisBoundaries::Periodic,
};
const WALLS: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Reflecting,
    upper: BoundaryCondition::Reflecting,
};

mod molecular;

fn mechanism(heat: f64, rate: f64, reference: f64) -> (MixtureGas, IdealGasKineticsV1) {
    let gas = make_gas(vec![
        source("A", R, 3.5, 0.0, reference),
        source("B", R, 3.5, 0.0, reference - heat),
    ]);
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("spatial-isomerization").unwrap(),
        vec![(SpeciesId::new("A").unwrap(), 1)],
        vec![(SpeciesId::new("B").unwrap(), 1)],
        ArrheniusRateV1::new(
            Temperature::new(1.0),
            ReactionProgressRateQuantityV1::new(rate),
            0.0,
            MolarEnergyQuantityV1::new(0.0),
        )
        .unwrap(),
    )
    .unwrap();
    let mechanism = IdealGasKineticsV1::new(
        gas.source().clone(),
        "manufactured-volumetric-reacting-gas",
        1,
        vec![reaction],
        1e-12,
    )
    .unwrap();
    (gas, mechanism)
}

fn chemistry() -> IsochoricKineticsSolve {
    IsochoricKineticsSolve {
        method: IsochoricKineticsMethod::Sdirk2,
        minimum_temperature_k: 0.001,
        maximum_temperature_k: 50.0,
        energy_tolerance_j: 1e-11,
        relative_tolerance: 2e-6,
        absolute_amount_tolerance_mol: 1e-13,
        entropy_tolerance_j_k: 1e-11,
        initial_step_s: 0.001,
        maximum_step_s: 1.0,
        maximum_steps: 10000,
        ..IsochoricKineticsSolve::default()
    }
}

fn field(
    gas: MixtureGas,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    initial: impl Fn([f64; 3]) -> ([f64; 2], [f64; 3], f64),
) -> MixtureEuler3d {
    let count: usize = grid.shape().iter().product();
    let mut states = vec![0.0; count * 6];
    for (i, q) in states.chunks_exact_mut(6).enumerate() {
        let (partial, velocity, pressure) = initial(grid.cell_center(i).unwrap());
        state_into(&gas, &partial, velocity, pressure, q).unwrap();
    }
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid).unwrap();
    MixtureEuler3d::new(gas, grid, states, boundaries, bytes, || false).unwrap()
}

fn coupled(flow: MixtureEuler3d, mechanism: IdealGasKineticsV1) -> ReactingEuler3d {
    let source = CellChemistry::new(flow.gas().clone(), mechanism, chemistry(), || false).unwrap();
    let bytes = ReactingEuler3d::required_workspace_bytes(&flow, &source).unwrap();
    ReactingEuler3d::with_chemistry(flow, source, bytes, || false).unwrap()
}

fn uniform() -> ReactingEuler3d {
    let (gas, mechanism) = mechanism(40.0, 500.0, 0.0);
    coupled(
        field(
            gas,
            CartesianGrid::new(2, 2, 2, 0.5, 0.5, 0.5).unwrap(),
            PERIODIC,
            |_| ([1.0, 0.0], [0.7, -0.4, 0.2], 1.0),
        ),
        mechanism,
    )
}

fn balance(report: &ReactiveStepReport3d) {
    close(report.mass_residual_kg, 0.0, 2e-11);
    close(report.absolute_energy_residual_j, 0.0, 2e-10);
    for value in report.conservation_residual.species_masses {
        close(value, 0.0, 2e-11);
    }
    for value in report.conservation_residual.momentum {
        close(value, 0.0, 2e-11);
    }
    close(report.conservation_residual.energy, 0.0, 2e-10);
    assert_eq!(report.elapsed, report.transport.elapsed);
}

#[test]
fn g1_heat_release_changes_wave_cfl_and_matches_original_homogeneous_reactor() {
    let mut flow = uniform();
    let volume = flow.flow().grid().cell_volume();
    let reference_source = flow.mechanism().clone();
    let report = flow
        .step(Time::new(0.2), ReactiveStepControl::default(), || false)
        .unwrap();
    balance(&report);
    assert!(
        report.coupling_retries > 0,
        "exothermic acoustic speed must shorten this step"
    );
    assert_eq!(report.chemistry.cell_solves, 16);
    assert_eq!(flow.flow().time(), report.elapsed);
    assert_eq!(flow.last_step(), Some(&report));
    assert!(report.chemistry.entropy_change_j_k > 0.0);
    let reference = advance_adiabatic(
        &reference_source,
        &[(SpeciesId::new("A").unwrap(), Amount::new(volume / R))],
        Volume::new(volume),
        Temperature::new(1.0),
        report.elapsed,
        chemistry(),
        || false,
    )
    .unwrap();
    for i in 0..flow.flow().cell_count() {
        let q = flow.flow().cell(i).unwrap();
        let p = q.primitive(flow.flow().gas(), || false).unwrap();
        close(p.temperature, reference.temperature().value(), 2e-5);
        close(
            q.partial_densities[1],
            reference.amounts()[1].value() * R / volume,
            2e-5,
        );
        for axis in 0..3 {
            close(p.velocity[axis], [0.7, -0.4, 0.2][axis], 2e-12);
        }
        assert!(p.temperature > 5.0);
    }
}

fn reacting_contact(n: usize) -> f64 {
    let (gas, mechanism) = mechanism(0.0, 2.0, 0.0);
    let dx = 1.0 / n as f64;
    let sinc = (std::f64::consts::PI * dx).sin() / (std::f64::consts::PI * dx);
    let average = sinc.powi(3);
    let velocity = [0.4, 0.25, -0.15];
    let exact = |x: [f64; 3], time: f64| {
        0.5 + 0.2
            * (-4.0 * time).exp()
            * average
            * (2.0 * std::f64::consts::PI * (x[0] + x[1] + x[2] - 0.5 * time)).sin()
    };
    let grid = CartesianGrid::new(n, n, n, dx, dx, dx).unwrap();
    let mut flow = coupled(
        field(gas, grid, PERIODIC, |x| {
            let a = exact(x, 0.0);
            ([a, 1.0 - a], velocity, 1.0)
        }),
        mechanism,
    );
    let end = 0.03;
    let mut steps = 0;
    while flow.flow().time().value() < end {
        let report = flow
            .step(
                Time::new(end - flow.flow().time().value()),
                ReactiveStepControl::default(),
                || false,
            )
            .unwrap();
        balance(&report);
        steps += 1;
        assert!(steps < 100);
    }
    let mut error = 0.0;
    for i in 0..flow.flow().cell_count() {
        let q = flow.flow().cell(i).unwrap();
        let p = q.primitive(flow.flow().gas(), || false).unwrap();
        close(p.pressure, 1.0, 2e-10);
        close(p.temperature, 1.0, 2e-10);
        error += (q.partial_densities[0] - exact(grid.cell_center(i).unwrap(), end)).abs();
    }
    error / flow.flow().cell_count() as f64
}

#[test]
fn g1_oblique_three_dimensional_advection_reaction_converges_to_exact_decay() {
    // Six points per phase put the MC extremum limiter outside its asymptotic
    // range: the oblique wave has only six distinct phase samples. Resolve the
    // wave before testing smooth spatial order against its exact cell averages.
    let coarse = reacting_contact(12);
    let fine = reacting_contact(24);
    let order = (coarse / fine).ln() / 2.0_f64.ln();
    eprintln!("3-D reacting contact L1: n12={coarse:e}, n24={fine:e}, order={order:e}");
    assert!(coarse < 0.01);
    assert!(
        order > 1.8,
        "coarse={coarse:e}, fine={fine:e}, order={order:e}"
    );
}

#[test]
fn g2_reacting_pocket_launches_compressive_waves_and_closes_chemical_energy_and_support_impulse() {
    let (gas, mechanism) = mechanism(24.0, 30.0, 0.0);
    let grid = CartesianGrid::new(6, 6, 6, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0).unwrap();
    let mut flow = coupled(
        field(
            gas,
            grid,
            Boundaries3d {
                x: WALLS,
                y: WALLS,
                z: WALLS,
            },
            |x| {
                let core = x.iter().map(|v| (v - 0.5).powi(2)).sum::<f64>() < 0.25 * 0.25;
                (if core { [1.0, 0.0] } else { [0.0, 1.0] }, [0.0; 3], 1.0)
            },
        ),
        mechanism,
    );
    let before = flow.flow().totals();
    let mut support = [0.0; 3];
    let end = 0.08;
    let mut steps = 0;
    while flow.flow().time().value() < end {
        let report = flow
            .step(
                Time::new(end - flow.flow().time().value()),
                ReactiveStepControl::default(),
                || false,
            )
            .unwrap();
        balance(&report);
        for side in BoundarySide::ALL {
            for impulse in flow.boundary_impulses(side) {
                assert_eq!(impulse.species_masses, [0.0, 0.0]);
                assert_eq!(impulse.energy, 0.0);
                for axis in 0..3 {
                    support[axis] += impulse.momentum[axis];
                }
            }
        }
        steps += 1;
        assert!(steps < 200);
    }
    let after = flow.flow().totals();
    close(after.mass(), before.mass(), 2e-12);
    close(
        after.absolute_energy(flow.flow().gas()).unwrap(),
        before.absolute_energy(flow.flow().gas()).unwrap(),
        2e-10,
    );
    for axis in 0..3 {
        close(
            after.momentum[axis] - before.momentum[axis] + support[axis],
            0.0,
            2e-12,
        );
    }
    assert!(after.species_masses[0] < before.species_masses[0]);
    let (mut peak_pressure, mut peak_speed, mut compression): (f64, f64, f64) = (0.0, 0.0, 0.0);
    for i in 0..flow.flow().cell_count() {
        let q = flow.flow().cell(i).unwrap();
        let p = q.primitive(flow.flow().gas(), || false).unwrap();
        peak_pressure = peak_pressure.max(p.pressure);
        peak_speed = peak_speed.max(p.velocity.iter().map(|v| v * v).sum::<f64>().sqrt());
        let outside = grid
            .cell_center(i)
            .unwrap()
            .iter()
            .map(|v| (v - 0.5).powi(2))
            .sum::<f64>()
            > 0.25 * 0.25;
        if outside {
            compression = compression.max(p.density - 1.0);
        }
    }
    eprintln!(
        "reacting pocket: steps={steps}, peak_p={peak_pressure:e}, speed={peak_speed:e}, outer_compression={compression:e}"
    );
    assert!(peak_pressure > 1.3);
    assert!(peak_speed > 0.05);
    assert!(compression > 0.01);
}

#[test]
fn g4_spatial_reacting_publication_is_atomic_and_all_field_memory_is_admitted() {
    let mut expected = uniform();
    let mut actual = uniform();
    for flow in [&mut expected, &mut actual] {
        flow.step(Time::new(0.003), ReactiveStepControl::default(), || false)
            .unwrap();
    }
    let before = actual.flow().states().to_vec();
    let before_time = actual.flow().time();
    let before_report = *actual.last_step().unwrap();
    let before_loads: Vec<_> = BoundarySide::ALL
        .into_iter()
        .flat_map(|side| {
            actual
                .boundary_impulses(side)
                .map(|v| (v.species_masses.to_vec(), v.momentum, v.energy))
        })
        .collect();
    let mut polls = 0;
    let reference = expected
        .step(Time::new(0.003), ReactiveStepControl::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut count = 0;
    let error = actual
        .step(Time::new(0.003), ReactiveStepControl::default(), || {
            count += 1;
            count >= polls
        })
        .unwrap_err();
    assert_eq!(error, ReactiveError::Cancelled);
    assert_eq!(actual.flow().states(), before);
    assert_eq!(actual.flow().time(), before_time);
    assert_eq!(actual.last_step(), Some(&before_report));
    let after_loads: Vec<_> = BoundarySide::ALL
        .into_iter()
        .flat_map(|side| {
            actual
                .boundary_impulses(side)
                .map(|v| (v.species_masses.to_vec(), v.momentum, v.energy))
        })
        .collect();
    assert_eq!(after_loads, before_loads);
    let replay = actual
        .step(Time::new(0.003), ReactiveStepControl::default(), || false)
        .unwrap();
    assert_eq!(replay, reference);
    assert_eq!(actual.flow().states(), expected.flow().states());

    let (gas, mechanism) = mechanism(0.0, 2.0, 0.0);
    let input = field(
        gas,
        CartesianGrid::new(2, 2, 2, 0.5, 0.5, 0.5).unwrap(),
        PERIODIC,
        |_| ([0.5, 0.5], [0.0; 3], 1.0),
    );
    let source = CellChemistry::new(input.gas().clone(), mechanism, chemistry(), || false).unwrap();
    let required = ReactingEuler3d::required_workspace_bytes(&input, &source).unwrap();
    assert_eq!(
        required,
        ReactingEuler3d::required_field_workspace_bytes(&input).unwrap() + source.workspace_bytes()
    );
    assert!(source.workspace_bytes() > 0);
    assert!(matches!(
        ReactingEuler3d::with_chemistry(input, source, required - 1, || false),
        Err(ReactiveError::Invalid(_))
    ));

    let mut refused = uniform();
    let original = refused.flow().states().to_vec();
    assert!(matches!(
        refused.step(
            Time::new(0.2),
            ReactiveStepControl {
                maximum_coupling_retries: 0,
                ..ReactiveStepControl::default()
            },
            || false
        ),
        Err(ReactiveError::RetryLimit { .. })
    ));
    assert_eq!(refused.flow().states(), original);
    assert_eq!(refused.flow().time().value(), 0.0);
    assert!(refused.last_step().is_none());
}
