use super::*;
use fs_evidence::ValidityDomain;
use fs_gas::dynamics::{Boundaries, Boundary, CaloricDomain, Euler1d};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Length, MolarMass, Pressure, Temperature};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K, required_nasa9_parameter_dims,
};

mod molecular;

pub(super) fn gas(cp_linear: f64) -> FrozenGas {
    // Manufactured NASA-9 source: R_specific=1, Cp/R=3.5+cp_linear*T.
    // No experimental material or constant-gamma claim for cp_linear != 0.
    let mut a = [0.0; 9];
    a[2] = 3.5;
    a[3] = cp_linear;
    let mut parameters = std::collections::BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(a) {
        parameters.insert(name.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 1.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()),
        law_version: 1,
        parameters,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 0.001, 50.0),
        sources: Vec::new(),
        provenance: Provenance {
            source: "Manufactured planar Euler calorics".into(),
            license: "test fixture".into(),
            artifact: None,
        },
    };
    let convention = StandardStateConventionV1::new(
        StandardStatePhaseV1::Gas,
        ReferenceEquationOfStateV1::IdealGas,
        Pressure::new(1.0),
        ElementalReferenceIdV1::new("planar-euler-fixture").unwrap(),
    )
    .unwrap();
    let species = Nasa9StandardStateModelV1::new(
        SpeciesId::new("test-gas").unwrap(),
        MolarMass::new(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K),
        convention,
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap();
    let source = FrozenIdealGasMixtureModelV1::new(
        vec![species],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
    )
    .unwrap();
    FrozenGas::new(
        source,
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap()
}

const PERIODIC: Boundaries2d = Boundaries2d {
    x: AxisBoundaries::Periodic,
    y: AxisBoundaries::Periodic,
};
const WALLS: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Reflecting,
    upper: BoundaryCondition::Reflecting,
};
const OPEN: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Transmissive,
    upper: BoundaryCondition::Transmissive,
};

fn field(
    gas: FrozenGas,
    grid: CartesianGrid,
    boundaries: Boundaries2d,
    primitive: impl Fn([f64; 2]) -> [f64; 4],
) -> Euler2d {
    let states = (0..grid.len())
        .map(|i| {
            let p = primitive(grid.cell_center(i).unwrap());
            State::from_primitive(&gas, p[0], [p[1], p[2]], p[3]).unwrap()
        })
        .collect();
    Euler2d::new(gas, grid, states, boundaries, || false).unwrap()
}

pub(super) fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn check_balance(report: &StepReport) {
    assert!(report.maximum_cfl <= 0.40000000000001);
    for residual in report.conservation_residual.array() {
        close(residual, 0.0, 4e-12);
    }
}

fn run(flow: &mut Euler2d, end: f64) {
    let mut steps = 0;
    while flow.time().value() < end {
        let report = flow
            .step(
                Time::new(end - flow.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        check_balance(&report);
        steps += 1;
        assert!(steps <= 2000, "bounded test duration");
    }
}

#[test]
fn g0_planar_periodic_uniform_state_and_physical_geometry() {
    let grid = CartesianGrid::new(7, 5, 0.12, 0.09, 0.4).unwrap();
    let mut flow = field(gas(0.4), grid, PERIODIC, |_| [1.1, 0.6, -0.3, 1.7]);
    let states = flow.states().to_vec();
    let initial = flow.totals();
    close(initial.mass, 1.1 * 7.0 * 5.0 * 0.12 * 0.09 * 0.4, 1e-14);
    close(grid.face_area(BoundarySide::XMin), 0.09 * 0.4, 1e-15);
    close(grid.face_area(BoundarySide::YMax), 0.12 * 0.4, 1e-15);
    run(&mut flow, 0.1);
    assert_eq!(flow.states(), states);
    assert_eq!(flow.totals(), initial);
    for q in flow.states() {
        let p = q.primitive(flow.gas(), || false).unwrap();
        close(p.pressure, 1.7, 1e-11);
        close(p.velocity[1], -0.3, 1e-14);
    }
}

#[test]
fn g1_oblique_contact_advects_with_second_order_refinement() {
    let end = 0.12;
    let mut errors = Vec::new();
    for n in [12, 24, 48] {
        let dx = 1.0 / n as f64;
        let z = std::f64::consts::PI * dx;
        let average = (z.sin() / z).powi(2);
        let exact = |x: [f64; 2], time: f64| {
            1.0 + 0.15 * average * (2.0 * std::f64::consts::PI * (x[0] + x[1] - 0.75 * time)).sin()
        };
        let grid = CartesianGrid::new(n, n, dx, dx, 0.6).unwrap();
        let mut flow = field(gas(0.0), grid, PERIODIC, |x| {
            [exact(x, 0.0), 0.5, 0.25, 1.0]
        });
        let before = flow.totals();
        run(&mut flow, end);
        let error = flow
            .states()
            .iter()
            .enumerate()
            .map(|(i, q)| {
                let p = q.primitive(flow.gas(), || false).unwrap();
                close(p.pressure, 1.0, 2e-10);
                (q.density - exact(grid.cell_center(i).unwrap(), end)).abs()
            })
            .sum::<f64>()
            / (n * n) as f64;
        errors.push(error);
        for (actual, expected) in flow.totals().array().into_iter().zip(before.array()) {
            close(actual, expected, 2e-12);
        }
    }
    eprintln!("oblique-contact density L1, 12/24/48 cells per axis: {errors:?}");
    assert!(errors[1] < 0.36 * errors[0]);
    assert!(errors[2] < 0.36 * errors[1]);
    assert!(errors[2] < 3e-4);
}

#[test]
fn g2_planar_sod_matches_one_dimensional_solution_and_axis_rotation() {
    let n = 48;
    let dx = 1.0 / n as f64;
    let g = gas(0.0);
    let initial = |position: f64| {
        if position < 0.5 {
            [1.0, 1.0]
        } else {
            [0.125, 0.1]
        }
    };
    let mut horizontal = field(
        g.clone(),
        CartesianGrid::new(n, 3, dx, 0.25, 0.7).unwrap(),
        Boundaries2d {
            x: OPEN,
            y: AxisBoundaries::Periodic,
        },
        |x| {
            let p = initial(x[0]);
            [p[0], 0.0, 0.0, p[1]]
        },
    );
    let mut vertical = field(
        g.clone(),
        CartesianGrid::new(3, n, 0.25, dx, 0.7).unwrap(),
        Boundaries2d {
            x: AxisBoundaries::Periodic,
            y: OPEN,
        },
        |x| {
            let p = initial(x[1]);
            [p[0], 0.0, 0.0, p[1]]
        },
    );
    let states = (0..n)
        .map(|i| {
            let p = initial((i as f64 + 0.5) * dx);
            g.conserved(p[0], 0.0, p[1]).unwrap()
        })
        .collect();
    let mut axial = Euler1d::new(
        g,
        Length::new(dx),
        0.75 * 0.7,
        states,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
        || false,
    )
    .unwrap();
    let end = 0.15;
    let mut steps = 0;
    while horizontal.time().value() < end {
        let x = horizontal
            .step(
                Time::new(end - horizontal.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        let y = vertical
            .step(x.elapsed, StepControl::default(), || false)
            .unwrap();
        let a = axial
            .step(x.elapsed, StepControl::default(), || false)
            .unwrap();
        close(y.elapsed.value(), x.elapsed.value(), 1e-14);
        close(a.elapsed.value(), x.elapsed.value(), 1e-14);
        check_balance(&x);
        check_balance(&y);
        steps += 1;
        assert!(steps < 1000);
    }
    for (i, a) in axial.states().iter().enumerate() {
        for j in 0..3 {
            let x = horizontal.states()[j * n + i];
            let y = vertical.states()[i * 3 + j];
            close(x.density, a.density, 2e-10);
            close(x.momentum[0], a.momentum, 2e-10);
            close(x.momentum[1], 0.0, 2e-12);
            close(x.energy, a.energy, 2e-10);
            close(y.density, x.density, 2e-10);
            close(y.momentum[1], x.momentum[0], 2e-10);
            close(y.momentum[0], 0.0, 2e-12);
            close(y.energy, x.energy, 2e-10);
        }
    }
    // The independent exact gamma-law Sod oracle is exercised by fs-gas's
    // canonical shock-tube tests; this cross-owner gate checks the complete
    // two-dimensional transport path against that same solver and source.
    let middle = horizontal.states()[n / 2]
        .primitive(horizontal.gas(), || false)
        .unwrap();
    assert!(middle.velocity[0] > 0.7 && middle.velocity[0] < 1.1);
    assert!(middle.pressure > 0.25 && middle.pressure < 0.4);
}

#[test]
fn g2_closed_wall_shock_transfers_equal_opposite_impulse_without_energy_loss() {
    let grid = CartesianGrid::new(20, 6, 0.05, 0.1, 0.4).unwrap();
    let mut flow = field(gas(0.4), grid, Boundaries2d { x: WALLS, y: WALLS }, |_| {
        [1.0, 2.0, 0.0, 1.0]
    });
    let initial = flow.totals();
    let mut impulse = [0.0; 2];
    let mut right_impulse = 0.0;
    let end = 0.04;
    while flow.time().value() < end {
        let report = flow
            .step(
                Time::new(end - flow.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        check_balance(&report);
        assert_eq!(report.boundary_change.mass, 0.0);
        assert_eq!(report.boundary_change.energy, 0.0);
        for side in SIDES {
            assert_eq!(flow.boundary_impulses(side).len(), grid.face_count(side));
            for face in flow.boundary_impulses(side) {
                assert_eq!(face.mass, 0.0);
                assert_eq!(face.energy, 0.0);
                impulse[0] += face.momentum[0];
                impulse[1] += face.momentum[1];
                if side == BoundarySide::XMax {
                    right_impulse += face.momentum[0];
                }
            }
        }
    }
    let final_state = flow.totals();
    close(final_state.mass, initial.mass, 2e-12);
    close(final_state.energy, initial.energy, 2e-12);
    for axis in 0..2 {
        close(
            final_state.momentum[axis] - initial.momentum[axis] + impulse[axis],
            0.0,
            2e-12,
        );
    }
    assert!(final_state.momentum[0] < initial.momentum[0]);
    assert!(
        right_impulse > 2.0 * end * 0.6 * 0.4,
        "reflected shock must raise wall pressure"
    );
    assert!(
        flow.states()[19]
            .primitive(flow.gas(), || false)
            .unwrap()
            .pressure
            > 2.0
    );
}

#[test]
fn g2_radial_pressure_blast_preserves_rotation_symmetry_and_closed_balances() {
    let n = 24;
    let grid = CartesianGrid::new(n, n, 1.0 / n as f64, 1.0 / n as f64, 0.5).unwrap();
    let mut flow = field(gas(0.0), grid, Boundaries2d { x: WALLS, y: WALLS }, |x| {
        let r2 = (x[0] - 0.5).powi(2) + (x[1] - 0.5).powi(2);
        [
            1.0,
            0.0,
            0.0,
            if r2 < 0.13_f64.powi(2) { 10.0 } else { 1.0 },
        ]
    });
    let before = flow.totals();
    run(&mut flow, 0.025);
    close(flow.totals().mass, before.mass, 3e-12);
    close(flow.totals().energy, before.energy, 3e-12);
    close(flow.totals().momentum[0], 0.0, 3e-12);
    close(flow.totals().momentum[1], 0.0, 3e-12);
    for j in 0..n {
        for i in 0..n {
            let q = flow.states()[j * n + i];
            let reflected = flow.states()[j * n + (n - i - 1)];
            let transpose = flow.states()[i * n + j];
            close(reflected.density, q.density, 3e-10);
            close(reflected.energy, q.energy, 3e-10);
            close(reflected.momentum[0], -q.momentum[0], 3e-10);
            close(reflected.momentum[1], q.momentum[1], 3e-10);
            close(transpose.density, q.density, 3e-10);
            close(transpose.energy, q.energy, 3e-10);
            close(transpose.momentum[0], q.momentum[1], 3e-10);
            close(transpose.momentum[1], q.momentum[0], 3e-10);
            let p = q.primitive(flow.gas(), || false).unwrap();
            assert!(p.density > 0.0 && p.pressure > 0.0 && p.temperature > 0.0);
        }
    }
    assert!(flow.states().iter().any(|q| q.density > 1.05));
    assert!(flow.states().iter().any(|q| q.momentum[0].abs() > 0.1));
}

fn accepted_wall_state() -> Euler2d {
    let mut flow = field(
        gas(0.4),
        CartesianGrid::new(4, 3, 0.25, 0.2, 0.5).unwrap(),
        Boundaries2d { x: WALLS, y: WALLS },
        |x| [1.0 + 0.1 * x[0], 0.2, -0.1, 1.0],
    );
    flow.step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
    flow
}

#[test]
fn g4_cancelled_and_failed_steps_preserve_accepted_state_and_wall_impulses() {
    let mut probe = accepted_wall_state();
    let mut total_polls = 0;
    probe
        .step(Time::new(0.001), StepControl::default(), || {
            total_polls += 1;
            false
        })
        .unwrap();
    for cutoff in [1, 64, total_polls / 2, total_polls - 1, total_polls] {
        let mut flow = accepted_wall_state();
        let cells = flow.states().to_vec();
        let time = flow.time();
        let report = *flow.last_step().unwrap();
        let impulses: Vec<_> = SIDES
            .into_iter()
            .map(|side| flow.boundary_impulses(side).to_vec())
            .collect();
        let mut polls = 0;
        assert_eq!(
            flow.step(Time::new(0.001), StepControl::default(), || {
                polls += 1;
                polls == cutoff
            }),
            Err(Error::Cancelled)
        );
        assert_eq!(flow.states(), cells);
        assert_eq!(flow.time(), time);
        assert_eq!(flow.last_step(), Some(&report));
        for (side, previous) in SIDES.into_iter().zip(&impulses) {
            assert_eq!(flow.boundary_impulses(side), previous);
        }
        assert!(
            flow.step(Time::new(0.0), StepControl::default(), || false)
                .is_err()
        );
        assert_eq!(flow.states(), cells);
        assert_eq!(flow.last_step(), Some(&report));
        // Cancellation only dirties private workspace; the next valid step
        // must reproduce the uninterrupted accepted result.
        flow.step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        assert_eq!(flow.states(), probe.states());
        assert_eq!(flow.last_step(), probe.last_step());
    }
    let mut refused = field(
        gas(0.0),
        CartesianGrid::new(2, 1, 0.5, 1.0, 1.0).unwrap(),
        PERIODIC,
        |x| {
            if x[0] < 0.5 {
                [1.0, -2.0, 30.0, 0.4]
            } else {
                [1.0, 2.0, -30.0, 0.4]
            }
        },
    );
    let initial = refused.states().to_vec();
    let controls = StepControl {
        reconstruction: Reconstruction::FirstOrder,
        maximum_halvings: 2,
        ..StepControl::default()
    };
    assert!(matches!(
        refused.step(Time::new(0.01), controls, || false),
        Err(Error::StepRetries { attempts: 3, .. })
    ));
    assert_eq!(refused.states(), initial);
    assert_eq!(refused.time().value(), 0.0);
    assert_eq!(refused.last_step(), None);
    for side in SIDES {
        assert!(
            refused
                .boundary_impulses(side)
                .iter()
                .all(|q| *q == Totals::default())
        );
    }
}

#[test]
fn g0_grid_and_state_admission_refuse_invalid_or_unrepresentable_extensives() {
    assert!(CartesianGrid::new(0, 2, 1.0, 1.0, 1.0).is_err());
    assert!(CartesianGrid::new(usize::MAX, 2, 1.0, 1.0, 1.0).is_err());
    assert!(CartesianGrid::new(2, 2, 1.0, f64::NAN, 1.0).is_err());
    let g = gas(0.0);
    let grid = CartesianGrid::new(1, 1, 1e100, 1e100, 1e100).unwrap();
    let q = State::from_primitive(&g, 1e10, [0.0, 0.0], 1e10).unwrap();
    assert!(
        Euler2d::new(g, grid, vec![q], PERIODIC, || false).is_err(),
        "public totals must stay finite after admission"
    );
}
