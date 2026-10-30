use super::*;
use crate::compressible::{AxisBoundaries, Euler2d, FrozenGas, SIDES, State};
use fs_evidence::ValidityDomain;
use fs_gas::multicomponent::{Boundaries, Boundary, CaloricDomain, MixtureEuler1d};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Length, MolarMass, Pressure, Temperature};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    ChargeVector, Composition, CompositionBasis, ElementId, ElementalMatrix,
    ElementalReferenceIdV1, Nasa9RegionV1, Nasa9StandardStateModelV1, ReferenceEquationOfStateV1,
    SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, required_nasa9_parameter_dims,
};

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

fn source(
    name: &str,
    mass: f64,
    cp: f64,
    linear: f64,
    formation: f64,
) -> Nasa9StandardStateModelV1 {
    // Manufactured NASA-9 laws, not measured material validation. Distinct
    // masses, heat capacities and references exercise the original source data.
    let coefficients = [0.0, 0.0, cp, linear, 0.0, 0.0, 0.0, formation, 0.0];
    let mut parameters = std::collections::BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
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
        sources: vec![],
        provenance: Provenance {
            source: "manufactured evolving planar mixture".into(),
            license: "fixture".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        SpeciesId::new(name).unwrap(),
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(1.0),
            ElementalReferenceIdV1::new("planar-mixture-fixture").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
fn gas(models: Vec<Nasa9StandardStateModelV1>) -> MixtureGas {
    let ids: Vec<_> = models.iter().map(|model| model.species().clone()).collect();
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1; ids.len()]],
        )
        .unwrap(),
        ChargeVector::new(ids.clone(), vec![0; ids.len()]).unwrap(),
        models,
    )
    .unwrap();
    MixtureGas::new(
        bank,
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap()
}
fn standard(linear: f64) -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.0, 30.0),
        source("B", 2.0 * R, 3.5, linear, -60.0),
        source("absent-C", 0.5 * R, 4.5, 0.2, 7.0),
    ])
}
fn field(
    gas: MixtureGas,
    grid: CartesianGrid,
    boundaries: Boundaries2d,
    initial: impl Fn([f64; 2]) -> ([f64; 3], [f64; 2], f64),
) -> MixtureEuler2d {
    let mut states = vec![0.0; grid.len() * 6];
    for (index, row) in states.chunks_exact_mut(6).enumerate() {
        let (densities, velocity, pressure) = initial(grid.cell_center(index).unwrap());
        state_into(&gas, &densities, velocity, pressure, row).unwrap();
    }
    let bytes = MixtureEuler2d::required_workspace_bytes(&gas, &grid).unwrap();
    MixtureEuler2d::new(gas, grid, states, boundaries, bytes, || false).unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}
fn balance(report: StepReport) {
    assert!(report.maximum_cfl <= 0.40000000000001);
    for residual in report.conservation_residual.species_masses {
        close(residual, 0.0, 3e-12);
    }
    for residual in report.conservation_residual.momentum {
        close(residual, 0.0, 3e-12);
    }
    close(report.conservation_residual.energy, 0.0, 3e-12);
    close(report.absolute_energy_residual_j, 0.0, 3e-10);
}
fn run(flow: &mut MixtureEuler2d, end: f64) {
    let mut steps = 0;
    while flow.time().value() < end {
        balance(
            flow.step(
                Time::new(end - flow.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap(),
        );
        steps += 1;
        assert!(steps < 500);
    }
}

#[test]
fn g1_oblique_species_contact_refines_and_updates_original_local_calorics() {
    let end = 0.1;
    let mut errors = Vec::new();
    for n in [12, 24] {
        let dx = 1.0 / n as f64;
        let z = std::f64::consts::PI * dx;
        let average = (z.sin() / z).powi(2);
        let exact = |x: [f64; 2], time: f64| {
            0.55 + 0.2 * average * (2.0 * std::f64::consts::PI * (x[0] + x[1] - 0.65 * time)).sin()
        };
        let grid = CartesianGrid::new(n, n, dx, dx, 0.6).unwrap();
        let mut flow = field(standard(0.0), grid, PERIODIC, |x| {
            let first = exact(x, 0.0);
            ([first, 1.0 - first, 0.0], [0.4, 0.25], 1.0)
        });
        let before = flow.totals();
        let initial = flow.states().to_vec();
        run(&mut flow, end);
        let mut error = 0.0;
        for index in 0..grid.len() {
            let q = flow.cell(index).unwrap();
            let p = q.primitive(flow.gas(), || false).unwrap();
            let a = q.partial_densities[0];
            let b = q.partial_densities[1];
            assert_eq!(q.partial_densities[2], 0.0);
            close(p.density, 1.0, 2e-11);
            close(p.pressure, 1.0, 2e-10);
            // M_A=R and M_B=2R: composition alone changes local T through
            // the original mixture EOS, without a user-supplied gamma/R.
            close(p.temperature, 1.0 / (a + 0.5 * b), 2e-10);
            close(p.sound_speed.powi(2), 1.4, 2e-10);
            error += (a - exact(grid.cell_center(index).unwrap(), end)).abs();
        }
        assert_ne!(flow.states(), initial);
        errors.push(error / grid.len() as f64);
        let after = flow.totals();
        for species in 0..3 {
            close(
                after.species_masses[species],
                before.species_masses[species],
                3e-12,
            );
        }
        close(
            after.absolute_energy(flow.gas()).unwrap(),
            before.absolute_energy(flow.gas()).unwrap(),
            3e-12,
        );
    }
    eprintln!("oblique species contact L1 on 12/24 grids: {errors:?}");
    assert!(errors[1] < 0.4 * errors[0]);
    assert!(errors[1] < 1.5e-3);
}

#[test]
fn g2_planar_species_shock_matches_original_axial_owner_and_rotation() {
    let n = 28;
    let dx = 1.0 / n as f64;
    let initial = |position: f64| {
        if position < 0.5 {
            ([0.9, 0.1, 0.0], 1.0)
        } else {
            ([0.025, 0.1, 0.0], 0.1)
        }
    };
    let g = standard(0.35);
    let mut horizontal = field(
        g.clone(),
        CartesianGrid::new(n, 3, dx, 0.2, 0.7).unwrap(),
        Boundaries2d {
            x: OPEN,
            y: AxisBoundaries::Periodic,
        },
        |x| {
            let (d, p) = initial(x[0]);
            (d, [0.0, 0.0], p)
        },
    );
    let mut vertical = field(
        g.clone(),
        CartesianGrid::new(3, n, 0.2, dx, 0.7).unwrap(),
        Boundaries2d {
            x: AxisBoundaries::Periodic,
            y: OPEN,
        },
        |x| {
            let (d, p) = initial(x[1]);
            (d, [0.0, 0.0], p)
        },
    );
    let states = (0..n)
        .map(|i| {
            let (d, p) = initial((i as f64 + 0.5) * dx);
            g.conserved(&d, 0.0, p).unwrap()
        })
        .collect();
    let mut axial = MixtureEuler1d::new(
        g,
        Length::new(dx),
        0.6 * 0.7,
        states,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
        || false,
    )
    .unwrap();
    while horizontal.time().value() < 0.06 {
        let report = horizontal
            .step(
                Time::new(0.06 - horizontal.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        balance(report);
        let rotated = vertical
            .step(report.elapsed, StepControl::default(), || false)
            .unwrap();
        balance(rotated);
        let one_dimensional = axial
            .step(report.elapsed, StepControl::default(), || false)
            .unwrap();
        close(rotated.elapsed.value(), report.elapsed.value(), 1e-14);
        close(
            one_dimensional.elapsed.value(),
            report.elapsed.value(),
            1e-14,
        );
    }
    for i in 0..n {
        let a = axial.cell(i).unwrap();
        for j in 0..3 {
            let x = horizontal.cell(j * n + i).unwrap();
            let y = vertical.cell(i * 3 + j).unwrap();
            for species in 0..3 {
                close(
                    x.partial_densities[species],
                    a.partial_densities[species],
                    3e-9,
                );
                close(
                    y.partial_densities[species],
                    x.partial_densities[species],
                    3e-10,
                );
            }
            close(x.momentum[0], a.momentum, 3e-9);
            close(x.momentum[1], 0.0, 1e-13);
            close(y.momentum[1], x.momentum[0], 3e-10);
            close(y.momentum[0], 0.0, 1e-13);
            close(x.energy, a.energy, 3e-9);
            close(y.energy, x.energy, 3e-10);
        }
    }
    let mut different_card = field(
        standard(1.5),
        *horizontal.grid(),
        horizontal.boundaries(),
        |x| {
            let (d, p) = initial(x[0]);
            (d, [0.0, 0.0], p)
        },
    );
    run(&mut different_card, 0.06);
    let middle = horizontal
        .cell(n / 2)
        .unwrap()
        .primitive(horizontal.gas(), || false)
        .unwrap();
    let changed = different_card
        .cell(n / 2)
        .unwrap()
        .primitive(different_card.gas(), || false)
        .unwrap();
    assert!(
        (middle.velocity[0] - changed.velocity[0]).abs() > 1e-3,
        "original NASA heat-capacity data must reach conserved shock evolution"
    );
}

#[test]
fn g0_same_single_source_preserves_frozen_cartesian_path() {
    let model = source("only", R, 3.5, 0.2, 1e17);
    let mixture = gas(vec![model.clone()]);
    let frozen = FrozenGas::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![model],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let grid = CartesianGrid::new(8, 6, 0.125, 1.0 / 6.0, 0.4).unwrap();
    let mut flat = vec![0.0; grid.len() * 4];
    let mut typed = Vec::new();
    for (index, row) in flat.chunks_exact_mut(4).enumerate() {
        let [x, y] = grid.cell_center(index).unwrap();
        let phase = 2.0 * std::f64::consts::PI;
        let density = 1.0 + 0.1 * (phase * x).sin();
        let velocity = [
            0.2 + 0.02 * (phase * y).sin(),
            -0.1 + 0.01 * (phase * x).cos(),
        ];
        let pressure = 1.0 + 0.1 * (phase * y).cos();
        state_into(&mixture, &[density], velocity, pressure, row).unwrap();
        typed.push(State::from_primitive(&frozen, density, velocity, pressure).unwrap());
    }
    let bytes = MixtureEuler2d::required_workspace_bytes(&mixture, &grid).unwrap();
    let mut evolving = MixtureEuler2d::new(mixture, grid, flat, PERIODIC, bytes, || false).unwrap();
    let mut fixed = Euler2d::new(frozen, grid, typed, PERIODIC, || false).unwrap();
    for _ in 0..3 {
        let a = evolving
            .step(Time::new(0.003), StepControl::default(), || false)
            .unwrap();
        let b = fixed
            .step(a.elapsed, StepControl::default(), || false)
            .unwrap();
        close(a.elapsed.value(), b.elapsed.value(), 1e-15);
    }
    for (index, q) in fixed.states().iter().enumerate() {
        let actual = evolving.cell(index).unwrap();
        close(actual.partial_densities[0], q.density, 3e-10);
        for axis in 0..2 {
            close(actual.momentum[axis], q.momentum[axis], 3e-10);
        }
        close(actual.energy, q.energy, 3e-10);
    }
}

#[test]
fn g2_closed_mixture_walls_report_full_vector_impulse_without_mass_or_heat_leak() {
    let grid = CartesianGrid::new(10, 7, 0.1, 0.12, 0.4).unwrap();
    let mut flow = field(
        standard(0.35),
        grid,
        Boundaries2d { x: WALLS, y: WALLS },
        |_| ([0.7, 0.3, 0.0], [1.2, -0.8], 1.0),
    );
    let initial = flow.totals();
    let mut outward = [0.0; 2];
    while flow.time().value() < 0.025 {
        let report = flow
            .step(
                Time::new(0.025 - flow.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        balance(report);
        assert_eq!(report.boundary_change.species_masses, [0.0; MAX_SPECIES]);
        assert_eq!(report.boundary_change.energy, 0.0);
        for side in SIDES {
            assert_eq!(flow.boundary_impulses(side).len(), grid.face_count(side));
            for face in flow.boundary_impulses(side) {
                assert_eq!(face.species_masses, &[0.0; 3]);
                assert_eq!(face.energy, 0.0);
                assert_eq!(face.momentum[1 - side.axis()], 0.0);
                for axis in 0..2 {
                    outward[axis] += face.momentum[axis];
                }
            }
        }
    }
    let after = flow.totals();
    for axis in 0..2 {
        close(
            after.momentum[axis] + outward[axis],
            initial.momentum[axis],
            2e-12,
        );
    }
    close(after.energy, initial.energy, 2e-12);
    close(
        after.absolute_energy(flow.gas()).unwrap(),
        initial.absolute_energy(flow.gas()).unwrap(),
        2e-12,
    );
    assert!(after.momentum[0] < initial.momentum[0]);
    assert!(after.momentum[1] > initial.momentum[1]);
}

#[test]
fn g4_total_component_budget_and_late_cancellation_are_atomic() {
    let grid = CartesianGrid::new(6, 5, 0.2, 0.15, 0.4).unwrap();
    let make = || {
        field(standard(0.2), grid, PERIODIC, |[x, y]| {
            let a = 0.45 + 0.1 * x.sin() * y.cos();
            ([a, 1.0 - a, 0.0], [0.3, -0.2], 1.0 + 0.02 * x.cos())
        })
    };
    let mut accepted = make();
    accepted
        .step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
    let bytes = accepted.workspace_bytes();
    let before = accepted.states().to_vec();
    assert!(matches!(
        MixtureEuler2d::new(
            accepted.gas().clone(),
            grid,
            before.clone(),
            PERIODIC,
            bytes - 1,
            || false
        ),
        Err(Error::Invalid(_))
    ));
    let mut excess_capacity = Vec::with_capacity(before.len() + 64);
    excess_capacity.extend_from_slice(&before);
    assert!(matches!(
        MixtureEuler2d::new(
            accepted.gas().clone(),
            grid,
            excess_capacity,
            PERIODIC,
            bytes,
            || false
        ),
        Err(Error::Invalid(_))
    ));
    assert!(matches!(
        MixtureEuler2d::new(
            accepted.gas().clone(),
            grid,
            before[..before.len() - 1].to_vec(),
            PERIODIC,
            bytes,
            || false
        ),
        Err(Error::Invalid(_))
    ));
    let mut probe = make();
    probe
        .step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
    let mut polls = 0;
    let expected = probe
        .step(Time::new(0.002), StepControl::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    let time = accepted.time();
    let receipt = *accepted.last_step().unwrap();
    let boundary = accepted.boundary_impulses.clone();
    let mut seen = 0;
    assert!(matches!(
        accepted.step(Time::new(0.002), StepControl::default(), || {
            seen += 1;
            seen == polls
        }),
        Err(Error::Cancelled)
    ));
    assert_eq!(seen, polls);
    assert_eq!(accepted.states(), before);
    assert_eq!(accepted.time(), time);
    assert_eq!(accepted.last_step(), Some(&receipt));
    assert_eq!(accepted.boundary_impulses, boundary);
    let replay = accepted
        .step(Time::new(0.002), StepControl::default(), || false)
        .unwrap();
    assert_eq!(replay, expected);
    assert_eq!(accepted.states(), probe.states());
    assert_eq!(accepted.boundary_impulses, probe.boundary_impulses);
}

#[test]
fn g4_source_refuses_unresolved_shear_heating_without_publishing_a_partial_step() {
    // Both initial thermal states lie inside the original NASA interval. The
    // HLL mixing state must include both kinetic components; its unresolved
    // counter-shear heat would exceed Tmax=50, regardless of time-step size.
    let grid = CartesianGrid::new(2, 2, 0.5, 0.5, 0.4).unwrap();
    let mut flow = field(standard(0.0), grid, PERIODIC, |[x, _]| {
        if x < 0.5 {
            ([1.0, 0.0, 0.0], [-2.0, 30.0], 0.4)
        } else {
            ([0.0, 1.0, 0.0], [2.0, -30.0], 0.4)
        }
    });
    let before = flow.states().to_vec();
    let boundary = flow.boundary_impulses.clone();
    let controls = StepControl {
        reconstruction: Reconstruction::FirstOrder,
        maximum_halvings: 2,
        ..StepControl::default()
    };
    let error = flow.step(Time::new(1e-4), controls, || false).unwrap_err();
    assert!(matches!(error, Error::StepRetries { attempts: 3, last }
        if matches!(*last, Error::Invalid("mixture HLL intermediate state outside source domain"))));
    assert_eq!(flow.states(), before);
    assert_eq!(flow.boundary_impulses, boundary);
    assert_eq!(flow.time().value(), 0.0);
    assert_eq!(flow.last_step(), None);
}
