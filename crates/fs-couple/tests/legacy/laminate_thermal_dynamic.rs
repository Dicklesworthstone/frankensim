//! G1/G3 transient thermoelastic checks on the original heat and MITC owners.
//! All cards are synthetic implementation oracles, not material measurements.

mod laminate_thermal_common;
#[path = "laminate_thermal_dynamic/viscosity.rs"]
mod viscosity;
use laminate_thermal_common::*;

use fs_couple::laminate_thermal::reciprocal::{
    LaminateThermoelasticControls, LaminateThermoelasticError, ReciprocalLayeredThermalShell,
    dynamic::DynamicLaminateState,
};
use fs_couple::laminate_thermal::{LaminateShellDiscretization, LayeredThermalShell};
use fs_material::laminate::ReferenceSurface;
use fs_material::thermal::transport::stack::transient::implicit::{
    LayeredThermalSolveControls, LayeredThermalWorkspace,
};
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::{ConstraintFrame, MitcBudget, NodalDirectors, ShellConstraint};

const T0: f64 = 300.0;
const ALPHA: f64 = 1e-5;

fn controls() -> LaminateThermoelasticControls {
    LaminateThermoelasticControls {
        maximum_iterations: 64,
        heat_relaxation: 1.0,
        cell_heat_tolerance_j: 1e-8,
        energy_tolerance_j: 1e-7,
    }
}

fn field_internal(state: &DynamicLaminateState) -> f64 {
    state
        .cell_potentials()
        .iter()
        .map(|cell| cell.internal_energy_j)
        .sum::<f64>()
        + state.transverse_shear_energy_j()
        + state.drilling_gauge_energy_j()
        + state.kinetic_energy_j()
}

/// The symmetric axial mode of this triangle is an exact scalar oscillator.
/// With u_x = q*x and all other nodal components fixed, k = E*V and
/// m = rho*integral(x^2 dV) = rho*V/2. Reflection in y exchanges the two
/// free nodes, so equal loads or velocities keep this mode invariant. Its
/// vorticity and transverse shear vanish; no drilling or bubble correction
/// enters either independent scalar coefficient.
fn axial_triangle() -> LaminateShellDiscretization {
    let mut d = discretization();
    d.mesh = ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, -0.6, 0.0], [1.0, 0.6, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    d.directors = NodalDirectors::new(&d.mesh, vec![[0.0, 0.0, 1.0]; 3], 0.9).unwrap();
    d.constraints = (0..3)
        .map(|node| ShellConstraint {
            node,
            mask: [node == 0, true, true, true, true, true],
            frame: ConstraintFrame::Global,
        })
        .collect();
    d.section_angles_rad = vec![0.0];
    d.budget = MitcBudget {
        max_nodes: 3,
        max_triangles: 1,
    };
    d
}

#[test]
fn g1_zero_expansion_preserves_native_inertia_and_the_scalar_midpoint_oscillator() {
    let cards = [card_with_coefficients(
        "axial inertia oracle",
        0.0,
        0.0,
        0.0,
        CP,
    )];
    let body = prepare_body(&cards, &[1], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let initial_thermal = initial(&bound, &[T0]);
    let shell = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        axial_triangle(),
        limits(),
        || false,
    )
    .unwrap();
    let equilibrium = ReciprocalLayeredThermalShell::try_new(
        &shell,
        &initial_thermal,
        &[0.0; 18],
        &[[280.0, 360.0]],
        0.02,
        controls(),
        || false,
    )
    .unwrap();
    let mut velocity = [0.0; 18];
    velocity[6] = 0.05;
    velocity[12] = 0.05;
    let mut coupled = equilibrium
        .begin_transient(&velocity, &[[0.0; 2]], 1e-7, || false)
        .unwrap();
    let volume = AREA * H;
    let mass = RHO * volume / 2.0;
    let stiffness = E * volume;
    let dt = 2e-4;
    let beta = dt * dt * stiffness / (4.0 * mass);
    let energy = 0.5 * mass * velocity[6] * velocity[6];
    close(coupled.state().kinetic_energy_j(), energy, 1e-12);
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    let mut thermal = initial_thermal;
    let (mut q, mut v) = (0.0, velocity[6]);
    for step_index in 1..=12 {
        let next_q = ((1.0 - beta) * q + dt * v) / (1.0 + beta);
        let next_v = 2.0 * (next_q - q) / dt - v;
        let step = coupled
            .advance_to(
                &thermal,
                step_index as f64 * dt,
                &[0.0; 18],
                conduction(),
                LayeredThermalSolveControls::default(),
                controls(),
                &mut workspace,
                || false,
            )
            .unwrap();
        for dof in [6, 12] {
            close(
                coupled.state().response().displacement_6n[dof],
                next_q,
                1e-11,
            );
            close(coupled.state().nodal_velocity_6n()[dof], next_v, 1e-8);
        }
        close(
            coupled.state().kinetic_energy_j(),
            0.5 * mass * next_v * next_v,
            1e-9,
        );
        close(
            coupled.state().cell_potentials()[0].free_energy_j,
            0.5 * stiffness * next_q * next_q,
            1e-9,
        );
        close(field_internal(coupled.state()), energy, 1e-8);
        close(coupled.state().temperatures_k()[0], T0, 1e-9);
        close(step.report.thermoelastic.reversible_heat_j, 0.0, 1e-8);
        close(
            step.report.thermoelastic.external_mechanical_work_j,
            0.0,
            1e-12,
        );
        close(
            step.report.thermoelastic.combined_energy_residual_j,
            0.0,
            1e-7,
        );
        close(step.thermal.mass_kg(), RHO * volume, 1e-12);
        thermal = step.thermal;
        (q, v) = (next_q, next_v);
    }
    assert!(q > 0.0);
    assert!(
        v < velocity[6],
        "the retained oscillator must actually slow down"
    );
}

#[test]
fn g1_adiabatic_inertial_traction_matches_the_independent_coupled_quadratic() {
    let cards = [card_with_coefficients(
        "inertial cooling oracle",
        0.0,
        ALPHA,
        0.0,
        CP,
    )];
    let body = prepare_body(&cards, &[1], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let thermal = initial(&bound, &[T0]);
    let original_cells = thermal.cells().clone();
    let shell = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        axial_triangle(),
        limits(),
        || false,
    )
    .unwrap();
    let equilibrium = ReciprocalLayeredThermalShell::try_new(
        &shell,
        &thermal,
        &[0.0; 18],
        &[[280.0, 360.0]],
        0.02,
        controls(),
        || false,
    )
    .unwrap();
    let mut coupled = equilibrium
        .begin_transient(&[0.0; 18], &[[0.0; 2]], 1e-7, || false)
        .unwrap();
    let dt = 1e-3;
    let force = 2400.0;
    let mut load = [0.0; 18];
    load[6] = 0.5 * force;
    load[12] = 0.5 * force;
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    let step = coupled
        .advance_to(
            &thermal,
            dt,
            &load,
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut workspace,
            || false,
        )
        .unwrap();
    let volume = AREA * H;
    let mass = RHO * volume / 2.0;
    let stiffness = E * volume;
    let capacity = RHO * volume * CP;
    let denominator = stiffness + 4.0 * mass / (dt * dt);
    // q=(F+k*alpha*dT)/D follows from midpoint momentum and qdot.
    // At nu=0 with transverse normal strain held at zero, the original
    // elastic internal energy is U=.5*k*q^2+k*alpha*T0*q+k*alpha^2*(T0^2-T^2).
    // C*dT=W-dU with W=.5*k*(q-alpha*dT)*q then gives this quadratic.
    let a = -stiffness * ALPHA * ALPHA + 0.5 * stiffness * stiffness * ALPHA * ALPHA / denominator;
    let b = capacity - 2.0 * stiffness * ALPHA * ALPHA * T0
        + stiffness * ALPHA * (0.5 * force + stiffness * ALPHA * T0) / denominator;
    let c = stiffness * ALPHA * force * T0 / denominator;
    let delta_t = -2.0 * c / (b + (b * b - 4.0 * a * c).sqrt());
    let temperature = T0 + delta_t;
    let q = (force + stiffness * ALPHA * delta_t) / denominator;
    let v = 2.0 * q / dt;
    let kinetic = 0.5 * mass * v * v;
    let internal = 0.5 * stiffness * q * q
        + stiffness * ALPHA * T0 * q
        + stiffness * ALPHA * ALPHA * (T0 * T0 - temperature * temperature);
    close(coupled.state().temperatures_k()[0], temperature, 3e-8);
    close(coupled.state().response().displacement_6n[6], q, 2e-11);
    close(coupled.state().response().displacement_6n[12], q, 2e-11);
    close(coupled.state().nodal_velocity_6n()[6], v, 4e-8);
    close(coupled.state().kinetic_energy_j(), kinetic, 1e-8);
    close(step.report.kinetic_energy_change_j, kinetic, 1e-8);
    close(
        coupled.state().cell_potentials()[0].internal_energy_j,
        internal,
        1e-7,
    );
    close(
        step.report.thermoelastic.external_mechanical_work_j,
        0.5 * force * q,
        1e-8,
    );
    close(step.report.thermoelastic.external_heat_j, 0.0, 1e-8);
    close(
        step.report.thermoelastic.combined_energy_residual_j,
        0.0,
        1e-7,
    );
    close(
        step.thermal.energy_j() - thermal.energy_j() + internal + kinetic,
        0.5 * force * q,
        1e-7,
    );
    assert!(temperature < T0);
    assert!(
        q < 0.2 * force / stiffness,
        "this short inertial step must not jump to static equilibrium"
    );
    assert!(step.report.thermoelastic.reversible_heat_j < 0.0);
    assert_eq!(thermal.cells(), &original_cells);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn g1_gradient_heating_and_internal_director_inertia_keep_one_atomic_energy_history() {
    let cards = [
        card_with_coefficients("dynamic lower laminate", 0.0, ALPHA, 1e-7, CP),
        card_with_coefficients("dynamic upper laminate", 0.0, 2.0 * ALPHA, 0.0, CP),
    ];
    let body = prepare_body(&cards, &[2, 1], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let thermal = initial_with_sources(&bound, &[300.0, 310.0, 340.0], &[0.0, 0.0, 200.0]);
    let original_cells = thermal.cells().clone();
    let shell = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        discretization(),
        limits(),
        || false,
    )
    .unwrap();
    let seed = || {
        ReciprocalLayeredThermalShell::try_new(
            &shell,
            &thermal,
            &[0.0; 24],
            &[[280.0, 360.0]; 3],
            0.02,
            controls(),
            || false,
        )
        .unwrap()
        .begin_transient(&[0.0; 24], &[[2.0, -1.0], [-0.5, 1.5]], 1e-7, || false)
        .unwrap()
    };
    let mut coupled = seed();
    let initial_displacement = coupled.state().response().displacement_6n.clone();
    let initial_velocity = coupled.state().nodal_velocity_6n().to_vec();
    let initial_field = field_internal(coupled.state());
    let initial_kinetic = coupled.state().kinetic_energy_j();
    // Integral_triangle (27*lambda0*lambda1*lambda2)^2 = area*81/280.
    // These initial velocities are purely the independently retained residual
    // director coordinates, so no nodal/cross-mass term contributes here.
    let rotary_inertia = RHO * H * H * H / 12.0;
    let expected_kinetic = 0.5
        * (AREA / 2.0)
        * rotary_inertia
        * (81.0 / 280.0)
        * (2.0_f64.powi(2) + 1.0 + 0.5_f64.powi(2) + 1.5_f64.powi(2));
    close(coupled.state().kinetic_energy_j(), expected_kinetic, 1e-11);
    let mut workspace = LayeredThermalWorkspace::new(3).unwrap();
    let target = 1e-4;
    assert!(matches!(
        coupled.advance_to(
            &thermal,
            target,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            LaminateThermoelasticControls {
                maximum_iterations: 1,
                ..controls()
            },
            &mut workspace,
            || false,
        ),
        Err(LaminateThermoelasticError::NotConverged { .. })
    ));
    close(coupled.state().time_s(), 0.0, 0.0);
    assert_eq!(
        coupled.state().response().displacement_6n,
        initial_displacement
    );
    assert_eq!(coupled.state().nodal_velocity_6n(), initial_velocity);
    close(coupled.state().kinetic_energy_j(), initial_kinetic, 0.0);

    let step = coupled
        .advance_to(
            &thermal,
            target,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(coupled.state().cell_potentials().len(), 3);
    assert!(
        (coupled.state().temperatures_k()[0] - coupled.state().temperatures_k()[2]).abs() > 1.0
    );
    assert!(
        step.report.bubble_material_work_j.abs() > 1e-9,
        "the retained internal director force must perform actual work"
    );
    assert!(
        coupled
            .state()
            .nodal_velocity_6n()
            .iter()
            .any(|v| v.abs() > 1e-6)
    );
    close(
        step.report.thermoelastic.external_heat_j,
        200.0 * target,
        1e-8,
    );
    close(
        step.report.thermoelastic.external_mechanical_work_j,
        0.0,
        1e-12,
    );
    close(
        step.report.thermoelastic.combined_energy_residual_j,
        0.0,
        1e-7,
    );
    close(
        step.thermal.energy_j() + field_internal(coupled.state())
            - thermal.energy_j()
            - initial_field,
        200.0 * target,
        2e-7,
    );
    assert_eq!(thermal.cells(), &original_cells);

    let accepted_displacement = coupled.state().response().displacement_6n.clone();
    let accepted_velocity = coupled.state().nodal_velocity_6n().to_vec();
    assert!(matches!(
        coupled.advance_to(
            &thermal,
            2.0 * target,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut workspace,
            || false,
        ),
        Err(LaminateThermoelasticError::Invalid(_))
    ));
    assert_eq!(
        coupled.state().response().displacement_6n,
        accepted_displacement
    );
    assert_eq!(coupled.state().nodal_velocity_6n(), accepted_velocity);
    close(coupled.state().time_s(), target, 0.0);

    // Count a successful proposal from a fresh owner, then cancel that same
    // deterministic proposal at its last checkpoint. The velocity/kinetic
    // history must remain accepted together with displacement and heat.
    let mut measured = seed();
    let mut measured_workspace = LayeredThermalWorkspace::new(3).unwrap();
    let mut checkpoints = 0;
    measured
        .advance_to(
            &thermal,
            target,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut measured_workspace,
            || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    let mut cancelled = seed();
    let mut cancelled_workspace = LayeredThermalWorkspace::new(3).unwrap();
    let mut checkpoint = 0;
    assert!(matches!(
        cancelled.advance_to(
            &thermal,
            target,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut cancelled_workspace,
            || {
                checkpoint += 1;
                checkpoint == checkpoints
            },
        ),
        Err(LaminateThermoelasticError::Cancelled)
    ));
    assert_eq!(checkpoint, checkpoints);
    close(cancelled.state().time_s(), 0.0, 0.0);
    assert_eq!(
        cancelled.state().response().displacement_6n,
        initial_displacement
    );
    assert_eq!(cancelled.state().nodal_velocity_6n(), initial_velocity);
    close(
        cancelled.state().kinetic_energy_j(),
        expected_kinetic,
        1e-11,
    );
    assert_eq!(thermal.cells(), &original_cells);
}
