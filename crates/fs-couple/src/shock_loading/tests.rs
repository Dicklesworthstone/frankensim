//! G0/G1/G3: source-owned reflected shocks, original pressure ports and modal motion.
#[cfg(feature = "experimental-compressible")]
mod molecular;
#[cfg(feature = "experimental-compressible")]
mod spatial;

use super::*;
use fs_evidence::ValidityDomain;
use fs_gas::dynamics::{CaloricDomain, FrozenGas, Reconstruction, StepControl};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_modal::time::{ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_qty::{Length, MolarMass, Pressure, Temperature, Time};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K, required_nasa9_parameter_dims,
};

fn gas() -> FrozenGas {
    // Manufactured SI source: R=1 J/(kg K), Cp/R=3.5; no empirical gas claim.
    let mut coefficients = [0.0; 9];
    coefficients[2] = 3.5;
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
        validity: ValidityDomain::unconstrained().with("T", 0.001, 100.0),
        sources: Vec::new(),
        provenance: Provenance {
            source: "manufactured shock loading fixture".into(),
            license: "test fixture".into(),
            artifact: None,
        },
    };
    let convention = StandardStateConventionV1::new(
        StandardStatePhaseV1::Gas,
        ReferenceEquationOfStateV1::IdealGas,
        Pressure::new(1.0),
        ElementalReferenceIdV1::new("shock-loading-fixture").unwrap(),
    )
    .unwrap();
    let source = Nasa9StandardStateModelV1::new(
        SpeciesId::new("manufactured-gas").unwrap(),
        MolarMass::new(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K),
        convention,
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap();
    FrozenGas::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![source],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(100.0)),
        || false,
    )
    .unwrap()
}

fn flow(velocity: f64, boundaries: Boundaries) -> Euler1d {
    let gas = gas();
    let states = vec![gas.conserved(1.0, velocity, 1.0).unwrap(); 40];
    Euler1d::new(gas, Length::new(0.025), 1.0, states, boundaries, || false).unwrap()
}
const WALLS: Boundaries = Boundaries::Ends {
    left: Boundary::Reflecting,
    right: Boundary::Reflecting,
};
fn controls() -> StepControl {
    StepControl {
        reconstruction: Reconstruction::FirstOrder,
        ..StepControl::default()
    }
}

fn mesh() -> ShellMesh {
    // One square at x=2 with normal +x, unit area and centroid (2,1/2,1/2).
    ShellMesh::new(
        vec![
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [2.0, 1.0, 1.0],
            [2.0, 0.0, 1.0],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
    .unwrap()
}
fn selection(side: ShellFace) -> Vec<OrientedShellFacet> {
    (0..2)
        .map(|triangle| OrientedShellFacet { triangle, side })
        .collect()
}
fn budget() -> SurfaceBudget {
    SurfaceBudget {
        geometry: ProfileBudget {
            max_nodes: 4,
            max_triangles: 2,
            max_feature_evaluations: 0,
        },
        maximum_patches: 2,
        maximum_preparation_terms: 100,
    }
}
fn surface(side: ShellFace) -> ShockLoadSurface {
    ShockLoadSurface::new(&mesh(), vec![selection(side)], budget(), || false).unwrap()
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-11 * a.abs().max(b.abs()).max(1.0),
        "{a} != {b}"
    );
}

#[test]
fn reflected_shock_uses_both_accepted_stages_and_closes_fluid_support_momentum() {
    let mut flow = flow(2.0, WALLS);
    let initial = flow.totals().momentum;
    let report = flow.step(Time::new(0.001), controls(), || false).unwrap();
    let left =
        WallPressureImpulse::from_euler1d(&flow, &report, End::Left, [1.0, 0.0, 0.0], 0.0).unwrap();
    let right = WallPressureImpulse::from_euler1d(&flow, &report, End::Right, [1.0, 0.0, 0.0], 0.0)
        .unwrap();
    // Exact constant-gamma one-sided Hugoniot solution at the initial trace.
    let gamma = 1.4_f64;
    let a = 2.0 / (gamma + 1.0);
    let b = (gamma - 1.0) / (gamma + 1.0);
    let velocity_squared = 4.0;
    let initial_pressure = (2.0 * a
        + velocity_squared
        + (velocity_squared * (4.0 * a * (1.0 + b) + velocity_squared)).sqrt())
        / (2.0 * a);
    assert!(right.mean_gas_pressure_pa() > 4.0);
    assert!(right.mean_gas_pressure_pa() < initial_pressure);
    let l = surface(ShellFace::Negative)
        .project(&[left], [0.0; 3], || false)
        .unwrap();
    let r = surface(ShellFace::Positive)
        .project(&[right], [0.0; 3], || false)
        .unwrap();
    close(
        flow.totals().momentum - initial + l.resultant_impulse_n_s[0] + r.resultant_impulse_n_s[0],
        0.0,
    );
    close(
        r.resultant_impulse_n_s[0],
        -report.right_boundary_change.momentum,
    );
    assert_eq!(r.resultant_impulse_n_s[1..], [0.0; 2]);
    close(r.angular_impulse_n_m_s[1], 0.5 * r.resultant_impulse_n_s[0]);
    close(
        r.angular_impulse_n_m_s[2],
        -0.5 * r.resultant_impulse_n_s[0],
    );
}

#[test]
fn exterior_subtraction_and_original_p1_virtual_work_are_exact() {
    let mut flow = flow(0.0, WALLS);
    let report = flow.step(Time::new(0.001), controls(), || false).unwrap();
    let impulse =
        WallPressureImpulse::from_euler1d(&flow, &report, End::Right, [1.0, 0.0, 0.0], 0.25)
            .unwrap();
    let load = surface(ShellFace::Positive)
        .project(&[impulse], [2.0, 0.5, 0.5], || false)
        .unwrap();
    close(load.resultant_impulse_n_s[0], 0.75 * report.elapsed.value());
    close(load.gas_impulse_n_s[0], report.elapsed.value());
    close(load.exterior_impulse_n_s[0], -0.25 * report.elapsed.value());
    for torque in load.angular_impulse_n_m_s {
        close(torque, 0.0);
    }
    let displacement: Vec<_> = mesh()
        .nodes
        .iter()
        .map(|p| [1.0 + 2.0 * p[1] + 3.0 * p[2], 0.0, 0.0])
        .collect();
    // Independent integral of the affine displacement over a unit square.
    close(
        load.virtual_work_j(&displacement).unwrap(),
        0.75 * (1.0 + 1.0 + 1.5),
    );
    let equilibrium =
        WallPressureImpulse::from_euler1d(&flow, &report, End::Right, [1.0, 0.0, 0.0], 1.0)
            .unwrap();
    let zero = surface(ShellFace::Positive)
        .project(&[equilibrium], [0.0; 3], || false)
        .unwrap();
    for value in zero.nodal_impulses_n_s.iter().flatten() {
        close(*value, 0.0);
    }
}

#[test]
fn original_modal_owner_receives_the_prescribed_shock_impulse_without_impulse_loss() {
    let mut flow = flow(1.0, WALLS);
    let report = flow.step(Time::new(0.001), controls(), || false).unwrap();
    let impulse =
        WallPressureImpulse::from_euler1d(&flow, &report, End::Right, [1.0, 0.0, 0.0], 1.0)
            .unwrap();
    let load = surface(ShellFace::Positive)
        .project(&[impulse], [0.0; 3], || false)
        .unwrap();
    let mass = 2.0_f64;
    let force: f64 = load
        .mean_nodal_forces_n
        .iter()
        .map(|f| f[0] / mass.sqrt())
        .sum();
    let mut body = ModalAcousticTimeModel::try_free_mass(
        48_000,
        mass,
        0.0,
        0.0,
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let frame = body.step_duration(&[force], load.duration_s).unwrap();
    let velocity = body.states()[0].velocity_m_sqrt_kg_per_s / mass.sqrt();
    close(mass * velocity, load.resultant_impulse_n_s[0]);
    close(frame.input_work_j, frame.total_modal_energy_j);
    assert!(velocity > 0.0);
    // This verifies the prescribed-load consumer only. The stationary Euler
    // wall is not moved; no coupled fluid/body energy closure is asserted.
}

#[test]
fn topology_geometry_boundary_clock_and_cancellation_refusals_are_atomic() {
    let mesh = mesh();
    assert!(
        ShockLoadSurface::new(
            &mesh,
            vec![
                selection(ShellFace::Positive),
                selection(ShellFace::Negative)
            ],
            budget(),
            || false
        )
        .is_err()
    );
    let mut short = budget();
    short.maximum_preparation_terms = 47;
    assert!(
        ShockLoadSurface::new(&mesh, vec![selection(ShellFace::Positive)], short, || false)
            .is_err()
    );
    assert_eq!(
        ShockLoadSurface::new(
            &mesh,
            vec![selection(ShellFace::Positive)],
            budget(),
            || true
        )
        .unwrap_err(),
        Error::Cancelled
    );
    let mut flow = flow(1.0, WALLS);
    let report = flow.step(Time::new(0.001), controls(), || false).unwrap();
    let impulse =
        WallPressureImpulse::from_euler1d(&flow, &report, End::Right, [1.0, 0.0, 0.0], 0.0)
            .unwrap();
    assert!(
        WallPressureImpulse::from_euler1d(&flow, &report, End::Right, [2.0, 0.0, 0.0], 0.0)
            .is_err()
    );
    assert!(
        WallPressureImpulse::from_euler1d(&flow, &report, End::Right, [1.0, 0.0, 0.0], -1.0)
            .is_err()
    );
    assert!(
        surface(ShellFace::Negative)
            .project(&[impulse], [0.0; 3], || false)
            .is_err()
    );
    let mut bigger = mesh.clone();
    for p in &mut bigger.nodes {
        p[1] *= 2.0;
    }
    let bigger = ShockLoadSurface::new(
        &bigger,
        vec![selection(ShellFace::Positive)],
        budget(),
        || false,
    )
    .unwrap();
    assert!(bigger.project(&[impulse], [0.0; 3], || false).is_err());
    assert_eq!(
        surface(ShellFace::Positive)
            .project(&[impulse], [0.0; 3], || true)
            .unwrap_err(),
        Error::Cancelled
    );
    let accepted = flow.states().to_vec();
    let time = flow.time();
    assert!(
        surface(ShellFace::Positive)
            .project(&[], [0.0; 3], || false)
            .is_err()
    );
    assert_eq!(flow.states(), accepted);
    assert_eq!(flow.time(), time);
    flow.step(Time::new(0.001), controls(), || false).unwrap();
    assert!(
        WallPressureImpulse::from_euler1d(&flow, &report, End::Right, [1.0, 0.0, 0.0], 0.0)
            .is_err()
    );
    let mut open = super::tests::flow(
        0.0,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
    );
    let report = open.step(Time::new(0.001), controls(), || false).unwrap();
    assert!(
        WallPressureImpulse::from_euler1d(&open, &report, End::Right, [1.0, 0.0, 0.0], 0.0)
            .is_err()
    );
}

#[test]
fn mixture_owner_retains_the_same_wall_pressure_and_forbids_transporting_ends() {
    use fs_gas::multicomponent::{MixtureEuler1d, MixtureGas};
    use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
    use fs_thermochem::{ChargeVector, ElementId, ElementalMatrix};
    let frozen = gas();
    let models = frozen.source().components().to_vec();
    let ids = models
        .iter()
        .map(|model| model.species().clone())
        .collect::<Vec<_>>();
    let source = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0]).unwrap(),
        models,
    )
    .unwrap();
    let mixture = MixtureGas::new(source, frozen.domain(), || false).unwrap();
    let states = vec![mixture.conserved(&[1.0], 1.0, 1.0).unwrap(); 40];
    let mut mixed =
        MixtureEuler1d::new(mixture, Length::new(0.025), 1.0, states, WALLS, || false).unwrap();
    let mut single = flow(1.0, WALLS);
    let a = single.step(Time::new(0.001), controls(), || false).unwrap();
    let b = mixed.step(Time::new(0.001), controls(), || false).unwrap();
    let frozen_load =
        WallPressureImpulse::from_euler1d(&single, &a, End::Right, [1.0, 0.0, 0.0], 1.0).unwrap();
    let mixed_load =
        WallPressureImpulse::from_mixture_euler1d(&mixed, &b, End::Right, [1.0, 0.0, 0.0], 1.0)
            .unwrap();
    close(
        frozen_load.gas_pressure_impulse_pa_s(),
        mixed_load.gas_pressure_impulse_pa_s(),
    );
    let load = surface(ShellFace::Positive)
        .project(&[mixed_load], [0.0; 3], || false)
        .unwrap();
    close(
        load.resultant_impulse_n_s[0],
        -b.right_boundary_change.momentum - b.elapsed.value(),
    );
    let mut corrupted = b;
    corrupted.right_boundary_change.species_masses[0] = 1e-8;
    assert!(
        WallPressureImpulse::from_mixture_euler1d(
            &mixed,
            &corrupted,
            End::Right,
            [1.0, 0.0, 0.0],
            0.0
        )
        .is_err()
    );
}

#[cfg(feature = "experimental-compressible")]
#[test]
fn nonuniform_cartesian_shock_footprint_preserves_its_structural_moment() {
    use fs_flux::compressible::{
        AxisBoundaries, Boundaries2d, BoundaryCondition, BoundarySide, CartesianGrid, Euler2d,
        State,
    };
    let gas = gas();
    let grid = CartesianGrid::new(4, 2, 0.5, 0.5, 1.0).unwrap();
    let states = (0..8)
        .map(|i| {
            State::from_primitive(&gas, 1.0, [1.0, 0.0], if i < 4 { 1.0 } else { 4.0 }).unwrap()
        })
        .collect();
    let reflecting = AxisBoundaries::Ends {
        lower: BoundaryCondition::Reflecting,
        upper: BoundaryCondition::Reflecting,
    };
    let mut flow = Euler2d::new(
        gas,
        grid,
        states,
        Boundaries2d {
            x: reflecting,
            y: reflecting,
        },
        || false,
    )
    .unwrap();
    let axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    assert!(WallPressureImpulse::from_euler2d(&flow, BoundarySide::XMax, 0, axes, 0.0).is_err());
    flow.step(Time::new(0.001), controls(), || false).unwrap();
    let impulses = (0..2)
        .map(|face| {
            WallPressureImpulse::from_euler2d(&flow, BoundarySide::XMax, face, axes, 0.0).unwrap()
        })
        .collect::<Vec<_>>();
    let mesh = ShellMesh::new(
        vec![
            [2.0, 0.0, 0.0],
            [2.0, 0.5, 0.0],
            [2.0, 1.0, 0.0],
            [2.0, 0.0, 1.0],
            [2.0, 0.5, 1.0],
            [2.0, 1.0, 1.0],
        ],
        vec![[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4]],
    )
    .unwrap();
    let patches = (0..2)
        .map(|face| {
            (2 * face..2 * face + 2)
                .map(|triangle| OrientedShellFacet {
                    triangle,
                    side: ShellFace::Positive,
                })
                .collect()
        })
        .collect();
    let budget = SurfaceBudget {
        geometry: ProfileBudget {
            max_nodes: 6,
            max_triangles: 4,
            max_feature_evaluations: 0,
        },
        maximum_patches: 2,
        maximum_preparation_terms: 120,
    };
    let surface = ShockLoadSurface::new(&mesh, patches, budget, || false).unwrap();
    let load = surface
        .project(&impulses, [2.0, 0.5, 0.5], || false)
        .unwrap();
    let raw = flow.boundary_impulses(BoundarySide::XMax);
    let [lower, upper] = [raw[0].momentum[0], raw[1].momentum[0]];
    assert!(upper > lower);
    close(load.resultant_impulse_n_s[0], lower + upper);
    // Independent moments at the two rectangular face centroids. Averaging
    // the wall pressure before projection would incorrectly erase this torque.
    close(load.angular_impulse_n_m_s[2], 0.25 * (lower - upper));
    assert!(load.angular_impulse_n_m_s[2] < 0.0);
    close(load.angular_impulse_n_m_s[1], 0.0);
    assert!(WallPressureImpulse::from_euler2d(&flow, BoundarySide::XMax, 2, axes, 0.0).is_err());
    assert!(
        WallPressureImpulse::from_euler2d(&flow, BoundarySide::XMax, 0, [axes[0]; 2], 0.0).is_err()
    );
}
