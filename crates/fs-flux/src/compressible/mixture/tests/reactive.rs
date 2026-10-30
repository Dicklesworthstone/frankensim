use super::*;
mod molecular;
use crate::compressible::mixture::reactive::{
    CellChemistry, ReactingEuler2d, ReactiveError, ReactiveStepControl, ReactiveStepReport,
};
use fs_gas::chemistry::{IsochoricKineticsMethod, IsochoricKineticsSolve, advance_adiabatic};
use fs_gas::multicomponent::reactive::{ReactingEuler1d, ReactiveStepControl as AxialControl};
use fs_qty::{Amount, Time, Volume};
use fs_thermochem::kinetics::{
    ArrheniusRateV1, ElementaryReactionV1, IdealGasKineticsV1, ReactionProgressRateQuantityV1,
};
use fs_thermochem::{MolarEnergyQuantityV1, ReactionId};

fn make_mechanism(heat: f64, rate: f64, product_cp: f64) -> IdealGasKineticsV1 {
    let source = gas(vec![
        source("A", R, 3.5, 0.0, 0.0),
        source("B", R, product_cp, 0.0, -heat),
        source("absent-C", R, 4.5, 0.0, 2.0),
    ]);
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("planar-isomerization").unwrap(),
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
    IdealGasKineticsV1::new(
        source.source().clone(),
        "manufactured-planar-reactions",
        1,
        vec![reaction],
        1e-12,
    )
    .unwrap()
}
fn chemistry() -> IsochoricKineticsSolve {
    IsochoricKineticsSolve {
        method: IsochoricKineticsMethod::Sdirk2,
        minimum_temperature_k: 0.001,
        maximum_temperature_k: 50.0,
        energy_tolerance_j: 1e-11,
        relative_tolerance: 1e-7,
        absolute_amount_tolerance_mol: 1e-14,
        entropy_tolerance_j_k: 1e-11,
        initial_step_s: 1e-4,
        maximum_step_s: 1.0,
        maximum_steps: 30000,
        ..IsochoricKineticsSolve::default()
    }
}
fn raw_field(
    mechanism: &IdealGasKineticsV1,
    grid: CartesianGrid,
    boundaries: Boundaries2d,
    initial: impl Fn([f64; 2]) -> ([f64; 3], [f64; 2], f64),
) -> MixtureEuler2d {
    field(
        MixtureGas::new(
            mechanism.thermodynamics().clone(),
            CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
            || false,
        )
        .unwrap(),
        grid,
        boundaries,
        initial,
    )
}
fn reactive(flow: MixtureEuler2d, mechanism: &IdealGasKineticsV1) -> ReactingEuler2d {
    let chemistry =
        CellChemistry::new(flow.gas().clone(), mechanism.clone(), chemistry(), || false).unwrap();
    let bytes = ReactingEuler2d::required_workspace_bytes(&flow, &chemistry).unwrap();
    let flow = ReactingEuler2d::new(flow, chemistry, bytes, || false).unwrap();
    assert_eq!(flow.workspace_bytes(), bytes);
    flow
}
fn uniform(mechanism: &IdealGasKineticsV1) -> ReactingEuler2d {
    reactive(
        raw_field(
            mechanism,
            CartesianGrid::new(4, 3, 0.25, 1.0 / 3.0, 0.7).unwrap(),
            PERIODIC,
            |_| ([1.0, 0.0, 0.0], [0.7, -0.3], 1.0),
        ),
        mechanism,
    )
}
fn coupled_balance(report: &ReactiveStepReport) {
    assert_eq!(report.elapsed, report.transport.elapsed);
    assert!(report.transport.maximum_cfl <= 0.40000000000001);
    for value in report.conservation_residual.species_masses {
        close(value, 0.0, 3e-11);
    }
    for value in report.conservation_residual.momentum {
        close(value, 0.0, 3e-11);
    }
    close(report.conservation_residual.energy, 0.0, 3e-10);
    close(report.mass_residual_kg, 0.0, 3e-11);
    close(report.absolute_energy_residual_j, 0.0, 3e-10);
}

#[test]
fn g1_exothermic_planar_gas_restarts_at_one_duration_and_retains_both_momenta() {
    let mechanism = make_mechanism(40.0, 500.0, 3.5);
    let mut flow = uniform(&mechanism);
    let report = flow
        .step(Time::new(0.2), ReactiveStepControl::default(), || false)
        .unwrap();
    assert!(
        report.coupling_retries > 0,
        "heat release must shorten the original CFL"
    );
    coupled_balance(&report);
    assert_eq!(flow.flow().time(), report.elapsed);
    assert_eq!(flow.last_step(), Some(&report));
    assert_eq!(report.chemistry.cell_solves, 2 * flow.flow().cell_count());
    assert!(report.chemistry.entropy_change_j_k > 0.0);
    assert!(report.chemistry.referenced_energy_change_j > 0.0);
    let volume = flow.flow().grid().cell_volume();
    let reference = advance_adiabatic(
        &mechanism,
        &[(SpeciesId::new("A").unwrap(), Amount::new(volume / R))],
        Volume::new(volume),
        Temperature::new(1.0),
        report.elapsed,
        chemistry(),
        || false,
    )
    .unwrap();
    for index in 0..flow.flow().cell_count() {
        let q = flow.flow().cell(index).unwrap();
        let p = q.primitive(flow.flow().gas(), || false).unwrap();
        assert_eq!(q.partial_densities[2], 0.0);
        close(p.velocity[0], 0.7, 2e-12);
        close(p.velocity[1], -0.3, 2e-12);
        close(p.temperature, reference.temperature().value(), 3e-6);
        close(
            q.partial_densities[1],
            reference.amounts()[1].value() * R / volume,
            3e-6,
        );
        assert!(p.temperature > 5.0);
    }

    // Only the original species' formation-energy DATA changes: it reaches
    // pressure and acoustics through this same solver and reaction mechanism.
    let mut hot = uniform(&mechanism);
    let mut neutral = uniform(&make_mechanism(0.0, 500.0, 3.5));
    for flow in [&mut hot, &mut neutral] {
        coupled_balance(
            &flow
                .step(Time::new(0.001), ReactiveStepControl::default(), || false)
                .unwrap(),
        );
    }
    let hot_p = hot
        .flow()
        .cell(0)
        .unwrap()
        .primitive(hot.flow().gas(), || false)
        .unwrap();
    let neutral_p = neutral
        .flow()
        .cell(0)
        .unwrap()
        .primitive(neutral.flow().gas(), || false)
        .unwrap();
    close(neutral_p.temperature, 1.0, 2e-10);
    assert!(hot_p.temperature > 2.0 * neutral_p.temperature);
    assert!(hot_p.sound_speed > neutral_p.sound_speed);
}

fn reacting_contact(n: usize) -> f64 {
    let mechanism = make_mechanism(0.0, 2.0, 3.5);
    let dx = 1.0 / n as f64;
    let z = std::f64::consts::PI * dx;
    let average = (z.sin() / z).powi(2);
    let exact = |point: [f64; 2], time: f64| {
        0.5 + 0.2
            * average
            * (-4.0 * time).exp()
            * (2.0 * std::f64::consts::PI * (point[0] + point[1] - 0.65 * time)).sin()
    };
    let grid = CartesianGrid::new(n, n, dx, dx, 0.6).unwrap();
    let mut flow = reactive(
        raw_field(&mechanism, grid, PERIODIC, |point| {
            let a = exact(point, 0.0);
            ([a, 1.0 - a, 0.0], [0.4, 0.25], 1.0)
        }),
        &mechanism,
    );
    let end = 0.08;
    let mut steps = 0;
    while flow.flow().time().value() < end {
        coupled_balance(
            &flow
                .step(
                    Time::new(end - flow.flow().time().value()),
                    ReactiveStepControl::default(),
                    || false,
                )
                .unwrap(),
        );
        steps += 1;
        assert!(steps < 100);
    }
    let mut error = 0.0;
    for index in 0..grid.len() {
        let q = flow.flow().cell(index).unwrap();
        let p = q.primitive(flow.flow().gas(), || false).unwrap();
        close(p.temperature, 1.0, 3e-10);
        close(p.pressure, 1.0, 3e-10);
        assert_eq!(q.partial_densities[2], 0.0);
        error += (q.partial_densities[0] - exact(grid.cell_center(index).unwrap(), end)).abs();
    }
    error / grid.len() as f64
}

#[test]
fn g1_oblique_reacting_contact_matches_analytic_decay_and_refines() {
    let coarse = reacting_contact(10);
    let fine = reacting_contact(20);
    eprintln!("oblique reacting contact L1 10/20 grids: {coarse:e} / {fine:e}");
    assert!(coarse < 4e-3);
    assert!(fine < 0.5 * coarse, "coarse={coarse:e},fine={fine:e}");
}

#[test]
fn g3_reacting_planar_shock_preserves_original_axial_owner_and_axis_rotation() {
    let mechanism = make_mechanism(6.0, 6.0, 4.5);
    let n = 18;
    let dx = 1.0 / n as f64;
    let initial = |x: f64| {
        if x < 0.5 {
            ([0.9, 0.1, 0.0], 1.0)
        } else {
            ([0.3, 0.2, 0.0], 0.3)
        }
    };
    let horizontal = raw_field(
        &mechanism,
        CartesianGrid::new(n, 2, dx, 0.2, 0.5).unwrap(),
        Boundaries2d {
            x: OPEN,
            y: AxisBoundaries::Periodic,
        },
        |point| {
            let (species, pressure) = initial(point[0]);
            (species, [0.15, 0.0], pressure)
        },
    );
    let gas = horizontal.gas().clone();
    let states = (0..n)
        .map(|i| {
            let (species, pressure) = initial((i as f64 + 0.5) * dx);
            gas.conserved(&species, 0.15, pressure).unwrap()
        })
        .collect();
    let axial = MixtureEuler1d::new(
        gas,
        Length::new(dx),
        0.1,
        states,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
        || false,
    )
    .unwrap();
    let mut axial = ReactingEuler1d::new(axial, mechanism.clone(), chemistry(), || false).unwrap();
    let mut horizontal = reactive(horizontal, &mechanism);
    let mut vertical = reactive(
        raw_field(
            &mechanism,
            CartesianGrid::new(2, n, 0.2, dx, 0.5).unwrap(),
            Boundaries2d {
                x: AxisBoundaries::Periodic,
                y: OPEN,
            },
            |point| {
                let (species, pressure) = initial(point[1]);
                (species, [0.0, 0.15], pressure)
            },
        ),
        &mechanism,
    );
    for _ in 0..3 {
        let dt = Time::new(0.001);
        let ar = axial.step(dt, AxialControl::default(), || false).unwrap();
        let hr = horizontal
            .step(dt, ReactiveStepControl::default(), || false)
            .unwrap();
        let vr = vertical
            .step(dt, ReactiveStepControl::default(), || false)
            .unwrap();
        assert_eq!(ar.elapsed, dt);
        assert_eq!(hr.elapsed, dt);
        assert_eq!(vr.elapsed, dt);
        coupled_balance(&hr);
        coupled_balance(&vr);
        for index in 0..n {
            let a = axial.flow().cell(index).unwrap();
            let h = horizontal.flow().cell(index).unwrap();
            let v = vertical.flow().cell(2 * index).unwrap();
            for species in 0..3 {
                close(
                    h.partial_densities[species],
                    a.partial_densities[species],
                    2e-10,
                );
                close(
                    v.partial_densities[species],
                    a.partial_densities[species],
                    2e-10,
                );
            }
            close(h.momentum[0], a.momentum, 2e-10);
            close(v.momentum[1], a.momentum, 2e-10);
            close(h.momentum[1], 0.0, 2e-10);
            close(v.momentum[0], 0.0, 2e-10);
            close(h.energy, a.energy, 2e-10);
            close(v.energy, a.energy, 2e-10);
        }
    }
}

fn wall_field(mechanism: &IdealGasKineticsV1) -> ReactingEuler2d {
    reactive(
        raw_field(
            mechanism,
            CartesianGrid::new(4, 3, 0.25, 1.0 / 3.0, 0.7).unwrap(),
            Boundaries2d { x: WALLS, y: WALLS },
            |point| {
                (
                    [0.8, 0.2, 0.0],
                    [0.1, -0.2],
                    0.8 + 0.4 * point[0] + 0.1 * point[1],
                )
            },
        ),
        mechanism,
    )
}
fn impulses(flow: &ReactingEuler2d) -> Vec<([f64; 3], [f64; 2], f64)> {
    SIDES
        .into_iter()
        .flat_map(|side| flow.boundary_impulses(side))
        .map(|q| {
            (
                [
                    q.species_masses[0],
                    q.species_masses[1],
                    q.species_masses[2],
                ],
                q.momentum,
                q.energy,
            )
        })
        .collect()
}

#[test]
fn g4_reacting_wall_impulse_and_final_cancel_preserve_the_complete_accepted_state() {
    let mechanism = make_mechanism(5.0, 8.0, 3.5);
    let mut reference = wall_field(&mechanism);
    reference
        .step(Time::new(0.001), ReactiveStepControl::default(), || false)
        .unwrap();
    let mut polls = 0;
    let expected = reference
        .step(Time::new(0.002), ReactiveStepControl::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    coupled_balance(&expected);
    assert!(polls > 100);
    let mut outward = [0.0; 2];
    for (species, momentum, energy) in impulses(&reference) {
        assert_eq!(species, [0.0; 3]);
        assert_eq!(energy, 0.0);
        for axis in 0..2 {
            outward[axis] += momentum[axis];
        }
    }
    for axis in 0..2 {
        close(
            expected.after.momentum[axis] - expected.before.momentum[axis] + outward[axis],
            0.0,
            3e-12,
        );
    }
    let mut flow = wall_field(&mechanism);
    flow.step(Time::new(0.001), ReactiveStepControl::default(), || false)
        .unwrap();
    let before = flow.flow().states().to_vec();
    let clock = flow.flow().time();
    let last = *flow.last_step().unwrap();
    let advection_last = *flow.flow().last_step().unwrap();
    let before_impulses = impulses(&flow);
    let mut attempted = 0;
    assert_eq!(
        flow.step(Time::new(0.002), ReactiveStepControl::default(), || {
            attempted += 1;
            attempted >= polls
        })
        .unwrap_err(),
        ReactiveError::Cancelled
    );
    assert_eq!(flow.flow().states(), before);
    assert_eq!(flow.flow().time(), clock);
    assert_eq!(flow.last_step(), Some(&last));
    assert_eq!(flow.flow().last_step(), Some(&advection_last));
    assert_eq!(impulses(&flow), before_impulses);
    let replay = flow
        .step(Time::new(0.002), ReactiveStepControl::default(), || false)
        .unwrap();
    assert_eq!(replay, expected);
    assert_eq!(flow.flow().states(), reference.flow().states());
    assert_eq!(impulses(&flow), impulses(&reference));
}

#[test]
fn g0_complete_workspace_source_and_chemistry_work_limits_are_admitted() {
    let mechanism = make_mechanism(10.0, 30.0, 3.5);
    let make_flow = || uniform(&mechanism).into_flow();
    let flow = make_flow();
    let prepared =
        CellChemistry::new(flow.gas().clone(), mechanism.clone(), chemistry(), || false).unwrap();
    let required = ReactingEuler2d::required_workspace_bytes(&flow, &prepared).unwrap();
    assert_eq!(
        required,
        flow.workspace_bytes()
            + MixtureEuler2d::required_workspace_bytes(flow.gas(), flow.grid()).unwrap()
            + prepared.workspace_bytes()
            + 2 * flow.boundary_impulses.len() * size_of::<f64>()
    );
    assert!(matches!(
        ReactingEuler2d::new(flow, prepared, required - 1, || false),
        Err(ReactiveError::Invalid(_))
    ));

    let flow = make_flow();
    let alternate_domain = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(0.002), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let alternate_controls = IsochoricKineticsSolve {
        minimum_temperature_k: 0.002,
        ..chemistry()
    };
    let other = CellChemistry::new(
        alternate_domain,
        mechanism.clone(),
        alternate_controls,
        || false,
    )
    .unwrap();
    assert!(matches!(
        ReactingEuler2d::new(flow, other, MAX_MIXTURE_WORKSPACE_BYTES, || false),
        Err(ReactiveError::Invalid(_))
    ));

    let mut limits = chemistry();
    limits.maximum_steps = 1;
    limits.initial_step_s = 1e-12;
    let flow = make_flow();
    let source = CellChemistry::new(flow.gas().clone(), mechanism, limits, || false).unwrap();
    let before = flow.states().to_vec();
    let bytes = ReactingEuler2d::required_workspace_bytes(&flow, &source).unwrap();
    let mut flow = ReactingEuler2d::new(flow, source, bytes, || false).unwrap();
    assert!(matches!(
        flow.step(
            Time::new(0.001),
            ReactiveStepControl {
                maximum_coupling_retries: 0,
                ..ReactiveStepControl::default()
            },
            || false
        ),
        Err(ReactiveError::RetryLimit { .. })
    ));
    assert_eq!(flow.flow().states(), before);
    assert_eq!(flow.flow().time().value(), 0.0);
    assert!(flow.last_step().is_none());
    assert!(flow.flow().last_step().is_none());
}
