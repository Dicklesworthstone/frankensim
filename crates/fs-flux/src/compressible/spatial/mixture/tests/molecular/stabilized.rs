use super::*;
use crate::compressible::spatial::mixture::molecular::{
    InertialWall3d, StabilizedDissipationControl, StabilizedMolecularReport, WallInertiaControl3d,
};
use fs_qty::Mass;

fn prepare_stable(
    flow: &MixtureEuler3d,
    law: ConstantTransport,
) -> MolecularTransport3d<ConstantTransport> {
    let source = SpatialTransport::new(flow.gas().clone(), law, || false).unwrap();
    let bytes = MolecularTransport3d::required_stabilized_workspace_bytes(flow, &source).unwrap();
    MolecularTransport3d::new(flow, source, bytes, || false).unwrap()
}

fn stable_balance(report: StabilizedMolecularReport) {
    assert!(report.maximum_scaled_parabolic_number <= 0.40000000000001);
    assert!(report.transport.maximum_mass_frame_residual_kg_m2_s < 1e-11);
    for value in report.transport.conservation_residual.species_masses {
        close(value, 0.0, 4e-11);
    }
    for value in report.transport.conservation_residual.momentum {
        close(value, 0.0, 4e-11);
    }
    close(report.transport.conservation_residual.energy, 0.0, 4e-11);
    close(report.transport.absolute_energy_residual_j, 0.0, 4e-9);
    assert!(report.transport.wall_inertia.is_none());
}

// Independently evaluate the published final Legendre stability polynomial,
// rather than replaying the production RKL stage coefficients.
fn legendre_amplification(stages: usize, argument: f64) -> f64 {
    let s = stages as f64;
    let b = (s * s + s - 2.0) / (2.0 * s * (s + 1.0));
    let x = 1.0 + 4.0 * argument / (s * s + s - 2.0);
    let (mut older, mut previous) = (1.0, x);
    for degree in 2..=stages {
        let j = degree as f64;
        let next = ((2.0 * j - 1.0) * x * previous - (j - 1.0) * older) / j;
        (older, previous) = (previous, next);
    }
    1.0 - b + b * previous
}

#[test]
fn g1_stabilized_three_axis_mode_matches_legendre_and_original_enthalpy_transport() {
    let n = 8;
    let h = 1.0 / n as f64;
    let diffusion = 0.03;
    let duration = 0.25;
    let grid = CartesianGrid::new(n, n, n, h, h, h).unwrap();
    let wave = |x: [f64; 3]| (2.0 * std::f64::consts::PI * (x[0] + x[1] + x[2])).sin();
    let mut flow = field(equal_mass_gas(), grid, PERIODIC, |x| {
        let a = 0.5 + 0.15 * wave(x);
        ([a, 1.0 - a, 0.0], [0.0; 3], 1.0)
    });
    let mut operator = prepare_stable(&flow, law(flow.gas(), diffusion, 0.0, 0.0, 0.0));
    let before_bytes = operator.workspace_bytes();
    let report = operator
        .advance_stabilized(
            &mut flow,
            Time::new(duration),
            StabilizedDissipationControl::default(),
            || false,
        )
        .unwrap();
    stable_balance(report);
    assert_eq!(report.transport.accepted_steps, 1);
    assert_eq!(report.transport.rejected_steps, 0);
    assert!(report.maximum_stages_used > 2);
    assert_eq!(report.operator_evaluations, report.maximum_stages_used + 1);
    assert!(report.transport.maximum_parabolic_number > 5.0 * 0.4);
    assert!(operator.workspace_bytes() > before_bytes);
    assert_eq!(flow.time().value(), 0.0);
    assert!(flow.last_step().is_none());
    let eigenvalue = 12.0 * diffusion * (std::f64::consts::PI * h).sin().powi(2) / h.powi(2);
    let amplification = legendre_amplification(report.maximum_stages_used, -duration * eigenvalue);
    assert!(amplification > 0.0 && amplification < 1.0);
    for index in 0..grid.len() {
        let q = flow.cell(index).unwrap();
        let p = q.primitive(flow.gas(), || false).unwrap();
        close(
            q.partial_densities[0],
            0.5 + 0.15 * amplification * wave(grid.cell_center(index).unwrap()),
            3e-11,
        );
        assert_eq!(q.partial_densities[2], 0.0);
        close(p.density, 1.0, 3e-12);
        // Different original Cv and formation offsets exercise the complete
        // species-referenced enthalpy flux, not a scalar-only diffusion path.
        close(p.temperature, 1.0, 3e-10);
        close(p.pressure, 1.0, 3e-10);
    }
    let retained_bytes = operator.workspace_bytes();
    let zero = operator
        .advance_stabilized(
            &mut flow,
            Time::new(0.0),
            StabilizedDissipationControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(zero.operator_evaluations, 0);
    assert_eq!(zero.maximum_stages_used, 0);
    assert_eq!(operator.workspace_bytes(), retained_bytes);
}

fn heated_channel() -> (MixtureEuler3d, MolecularTransport3d<ConstantTransport>) {
    let flow = field(
        equal_mass_gas(),
        CartesianGrid::new(3, 4, 2, 1.0 / 3.0, 0.25, 0.5).unwrap(),
        Boundaries3d {
            y: WALLS,
            ..PERIODIC
        },
        |_| ([1.0, 0.0, 0.0], [0.05, 0.0, -0.03], 1.0),
    );
    let operator = prepare_stable(&flow, law(flow.gas(), 0.002, 0.03, 0.01, 0.04))
        .with_wall(
            BoundarySide::YMax,
            MolecularWall3d::NoSlip {
                velocity: [0.2, 0.0, 0.1],
                thermal: WallThermalCondition::Isothermal(Temperature::new(2.0)),
            },
        )
        .unwrap();
    (flow, operator)
}

#[test]
fn g1_stabilized_static_heated_wall_closes_signed_rkl_boundary_transfers() {
    let (mut flow, mut operator) = heated_channel();
    let report = operator
        .advance_stabilized(
            &mut flow,
            Time::new(0.2),
            StabilizedDissipationControl::default(),
            || false,
        )
        .unwrap();
    stable_balance(report);
    assert!(report.maximum_stages_used > 2);
    assert!(report.transport.boundary_change.energy > 0.0);
    assert!(report.transport.boundary_change.momentum[0] > 0.0);
    assert!(report.transport.boundary_change.momentum[2] > 0.0);
    for impulse in operator.boundary_impulses(BoundarySide::YMax) {
        assert_eq!(impulse.species_masses, [0.0; 3]);
        assert!(impulse.momentum[0] < 0.0 && impulse.momentum[2] < 0.0);
        assert!(impulse.energy < 0.0);
    }
    let outward: f64 = BoundarySide::ALL
        .into_iter()
        .flat_map(|side| operator.boundary_impulses(side))
        .map(|q| q.energy)
        .sum();
    close(
        report.transport.after.energy - report.transport.before.energy + outward,
        0.0,
        4e-11,
    );
    let published = *operator.last_step().unwrap();
    let loads = impulses(&operator);
    let state = flow.states().to_vec();
    let limited = StabilizedDissipationControl {
        maximum_operator_evaluations: 1,
        ..StabilizedDissipationControl::default()
    };
    assert!(matches!(
        operator.advance_stabilized(&mut flow, Time::new(0.2), limited, || false),
        Err(DissipationError::OperatorWorkLimit {
            operator_evaluations: 1,
            ..
        })
    ));
    assert_eq!(flow.states(), state);
    assert_eq!(*operator.last_step().unwrap(), published);
    assert_eq!(impulses(&operator), loads);
}

#[test]
fn g4_stabilized_late_cancellation_and_lazy_storage_refuse_without_publication() {
    let (mut reference, mut reference_operator) = heated_channel();
    let mut polls = 0;
    reference_operator
        .advance_stabilized(
            &mut reference,
            Time::new(0.2),
            StabilizedDissipationControl::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 10);
    let (mut flow, mut operator) = heated_channel();
    let initial = flow.states().to_vec();
    let initial_loads = impulses(&operator);
    let mut attempted = 0;
    assert!(matches!(
        operator.advance_stabilized(
            &mut flow,
            Time::new(0.2),
            StabilizedDissipationControl::default(),
            || {
                attempted += 1;
                attempted >= polls - 2
            },
        ),
        Err(DissipationError::Cancelled)
    ));
    assert_eq!(flow.states(), initial);
    assert_eq!(flow.time().value(), 0.0);
    assert!(operator.last_step().is_none());
    assert_eq!(impulses(&operator), initial_loads);
    let source = SpatialTransport::new(
        flow.gas().clone(),
        law(flow.gas(), 0.002, 0.03, 0.01, 0.04),
        || false,
    )
    .unwrap();
    let required =
        MolecularTransport3d::required_stabilized_workspace_bytes(&flow, &source).unwrap();
    let mut bounded = MolecularTransport3d::new(&flow, source, required - 1, || false).unwrap();
    let bytes = bounded.workspace_bytes();
    assert!(matches!(
        bounded.advance_stabilized(
            &mut flow,
            Time::new(0.1),
            StabilizedDissipationControl::default(),
            || false,
        ),
        Err(DissipationError::Invalid(
            "volumetric molecular total workspace budget"
        ))
    ));
    assert_eq!(bounded.workspace_bytes(), bytes);
    assert_eq!(flow.states(), initial);
    assert!(bounded.last_step().is_none());
}

#[test]
fn g0_stabilized_finite_wall_inertia_requires_its_original_reciprocal_solver() {
    let (mut flow, _) = heated_channel();
    let source = SpatialTransport::new(
        flow.gas().clone(),
        law(flow.gas(), 0.002, 0.03, 0.01, 0.04),
        || false,
    )
    .unwrap();
    let bytes = MolecularTransport3d::required_stabilized_workspace_bytes(&flow, &source).unwrap()
        + WallInertiaControl3d::solver_workspace_bytes();
    let mut operator = MolecularTransport3d::new(&flow, source, bytes, || false)
        .unwrap()
        .with_wall(
            BoundarySide::YMax,
            MolecularWall3d::NoSlip {
                velocity: [0.2, 0.0, 0.1],
                thermal: WallThermalCondition::Adiabatic,
            },
        )
        .unwrap()
        .with_wall_inertia(
            &[InertialWall3d {
                side: BoundarySide::YMax,
                mass: Mass::new(1.0),
            }],
            WallInertiaControl3d::default(),
        )
        .unwrap();
    let initial = flow.states().to_vec();
    assert!(matches!(
        operator.advance_stabilized(
            &mut flow,
            Time::new(0.1),
            StabilizedDissipationControl::default(),
            || false,
        ),
        Err(DissipationError::Invalid(
            "stabilized molecular transport requires prescribed walls without finite inertia"
        ))
    ));
    assert_eq!(flow.states(), initial);
    assert!(operator.last_step().is_none());
}
