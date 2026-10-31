//! Reciprocal heat/plane-stress shell checks on original physical owners.
//! Synthetic source cards; analytical implementation oracles, not measurements.

mod laminate_thermal_common;
use laminate_thermal_common::*;

use fs_couple::laminate_thermal::LayeredThermalShell;
use fs_couple::laminate_thermal::reciprocal::{
    LaminateThermoelasticControls, LaminateThermoelasticError, ReciprocalLaminateState,
    ReciprocalLayeredThermalShell,
};
use fs_material::laminate::ReferenceSurface;
use fs_material::thermal::transport::stack::transient::implicit::{
    LayeredThermalSolveControls, LayeredThermalWorkspace,
};
use fs_plate::shell::mitc::{ConstraintFrame, ShellConstraint};

const ALPHA: f64 = 1e-5;
const T0: f64 = 300.0;

fn controls() -> LaminateThermoelasticControls {
    LaminateThermoelasticControls {
        maximum_iterations: 64,
        heat_relaxation: 1.0,
        cell_heat_tolerance_j: 1e-8,
        energy_tolerance_j: 1e-7,
    }
}
fn field_internal(state: &ReciprocalLaminateState) -> f64 {
    state
        .cell_potentials()
        .iter()
        .map(|p| p.internal_energy_j)
        .sum::<f64>()
        + state.transverse_shear_energy_j()
        + state.drilling_gauge_energy_j()
}

#[test]
fn g1_adiabatic_axial_traction_cools_by_the_independent_finite_step_formula() {
    let cards = [card_with_coefficients(
        "adiabatic sheet",
        0.0,
        ALPHA,
        0.0,
        CP,
    )];
    let body = prepare_body(&cards, &[1], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let thermal = initial(&bound, &[T0]);
    let shell = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        discretization(),
        limits(),
        || false,
    )
    .unwrap();
    let mut coupled = ReciprocalLayeredThermalShell::try_new(
        &shell,
        &thermal,
        &[0.0; 24],
        &[[280.0, 360.0]],
        0.02,
        controls(),
        || false,
    )
    .unwrap();
    let sigma = 2e5;
    // Exact consistent tractions on the two vertical edges. Corner 0 remains
    // the original rigid-motion support; its applied force does zero work.
    let nodal_force = sigma * H * 0.6 / 2.0;
    let mut load = [0.0; 24];
    for (node, sign) in [(0, -1.0), (1, 1.0), (2, -1.0), (3, 1.0)] {
        load[6 * node] = sign * nodal_force;
    }
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    let step = coupled
        .advance_to(
            &thermal,
            1.0,
            &load,
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut workspace,
            || false,
        )
        .unwrap();
    let volume = AREA * H;
    let capacity = volume * RHO * CP;
    let coupling = 0.5 * volume * ALPHA * sigma;
    // C dT = -.5 V alpha sigma (T1+T0), from endpoint-average
    // mechanical work and U=Psi+T*S at the same compatible strain.
    let expected = T0 * (capacity - coupling) / (capacity + coupling);
    close(coupled.state().shell().temperatures_k()[0], expected, 3e-8);
    assert!(expected < T0);
    close(
        coupled.state().shell().response().displacement_6n[6],
        sigma / E + ALPHA * (expected - T0),
        2e-11,
    );
    close(
        coupled.state().shell().response().displacement_6n[13],
        0.6 * ALPHA * (expected - T0),
        2e-11,
    );
    close(step.report.external_heat_j, 0.0, 1e-8);
    assert!(step.report.reversible_heat_j < 0.0);
    close(
        step.thermal.energy_j() + field_internal(coupled.state()) - thermal.energy_j(),
        step.report.external_mechanical_work_j,
        1e-7,
    );
    close(step.report.combined_energy_residual_j, 0.0, 1e-7);
    close(step.report.mechanical_work_residual_j, 0.0, 1e-7);
    close(step.thermal.mass_kg(), RHO * AREA * H, 1e-12);
    assert_eq!(thermal.time_s(), 0.0);
    close(
        bound.cell_temperatures(&thermal, || Ok(())).unwrap()[0],
        T0,
        0.0,
    );
}

#[test]
fn g1_restrained_heating_uses_caloric_plus_elastic_capacity_with_zero_mechanical_work() {
    let cards = [card_with_coefficients(
        "restrained sheet",
        0.0,
        ALPHA,
        0.0,
        CP,
    )];
    let body = prepare_body(&cards, &[1], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let input_heat = 6000.0;
    let thermal = initial_with_sources(&bound, &[T0], &[input_heat]);
    let mut mesh = discretization();
    // Restrain both in-plane strains everywhere, retaining native bending
    // unknowns and the original corner support for a nonsingular solve.
    for node in 1..4 {
        mesh.constraints.push(ShellConstraint {
            node,
            mask: [true, true, false, false, false, true],
            frame: ConstraintFrame::Global,
        });
    }
    let shell = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        mesh,
        limits(),
        || false,
    )
    .unwrap();
    let mut coupled = ReciprocalLayeredThermalShell::try_new(
        &shell,
        &thermal,
        &[0.0; 24],
        &[[280.0, 360.0]],
        0.02,
        controls(),
        || false,
    )
    .unwrap();
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    let step = coupled
        .advance_to(
            &thermal,
            1.0,
            &[0.0; 24],
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut workspace,
            || false,
        )
        .unwrap();
    let volume = AREA * H;
    let c = volume * RHO * CP;
    let b = volume * E * ALPHA * ALPHA;
    // Both normal directions are restrained, nu=0: U_el=B(T0²-T²).
    // Thus (C-2BT0)dT-B*dT²=Q. Use the stable quadratic root.
    let effective = c - 2.0 * b * T0;
    let delta =
        2.0 * input_heat / (effective + (effective * effective - 4.0 * b * input_heat).sqrt());
    let expected = T0 + delta;
    let actual = coupled.state().shell().temperatures_k()[0];
    close(actual, expected, 3e-8);
    assert!(
        actual > T0 + input_heat / c,
        "elastic entropy changes the effective constrained capacity"
    );
    close(
        coupled.state().cell_potentials()[0].internal_energy_j,
        b * (T0 * T0 - actual * actual),
        1e-7,
    );
    close(
        coupled.state().cell_potentials()[0].heat_capacity_j_k,
        -2.0 * b * actual,
        1e-10,
    );
    close(step.report.external_mechanical_work_j, 0.0, 1e-12);
    close(step.report.external_heat_j, input_heat, 1e-7);
    close(step.report.combined_energy_residual_j, 0.0, 1e-7);
    assert!(step.report.reversible_heat_j > 0.0);
    // The accepted reciprocal temperature/stress reaches the public native
    // stability path. For unit corner w on this 1 x .6 two-triangle patch,
    // integral |grad(w)|^2 dA = .3*(1 + 1/.6^2).
    let prestress = coupled
        .state()
        .shell()
        .membrane_prestress(
            fs_plate::shell::mitc::thermal::prestress::MembranePrestressLimits {
                maximum_displacement_gradient: 0.05,
                maximum_rotation_rad: 0.05,
            },
            || false,
        )
        .unwrap();
    let mut corner = [0.0; 24];
    corner[20] = 1.0;
    let direction = prestress.project_load(&corner).unwrap();
    let geometric_work: f64 = direction
        .iter()
        .enumerate()
        .map(|(row, &u)| {
            let (columns, values) = prestress.geometric_stiffness().row(row);
            u * columns
                .iter()
                .zip(values)
                .map(|(&col, &value)| value * direction[col])
                .sum::<f64>()
        })
        .sum();
    let expected_geometric = -E * H * ALPHA * (actual - T0) * 0.3 * (1.0 + 1.0 / (0.6 * 0.6));
    close(geometric_work, expected_geometric, 1e-7);
    assert!(geometric_work < 0.0);
    for node in 0..4 {
        close(
            coupled.state().shell().response().displacement_6n[6 * node],
            0.0,
            0.0,
        );
        close(
            coupled.state().shell().response().displacement_6n[6 * node + 1],
            0.0,
            0.0,
        );
    }
}

#[test]
fn g1_spatial_bending_and_native_transverse_work_close_energy_and_preserve_atomic_history() {
    let cards = [
        card_with_coefficients("gradient lower", 0.0, ALPHA, 1e-7, CP),
        card_with_coefficients("gradient upper", 0.0, 2.0 * ALPHA, 0.0, CP),
    ];
    let body = prepare_body(&cards, &[2, 1], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let thermal = initial(&bound, &[300.0, 310.0, 340.0]);
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
    let ranges = [[280.0, 360.0]; 3];
    let mut coupled = ReciprocalLayeredThermalShell::try_new(
        &shell,
        &thermal,
        &[0.0; 24],
        &ranges,
        0.02,
        controls(),
        || false,
    )
    .unwrap();
    let initial_field = field_internal(coupled.state());
    let initial_displacement = coupled.state().shell().response().displacement_6n.clone();
    let mut load = [0.0; 24];
    load[20] = -1.0;
    let mut workspace = LayeredThermalWorkspace::new(3).unwrap();
    assert!(matches!(
        coupled.advance_to(
            &thermal,
            0.5,
            &load,
            conduction(),
            LayeredThermalSolveControls::default(),
            LaminateThermoelasticControls {
                maximum_iterations: 1,
                ..controls()
            },
            &mut workspace,
            || false
        ),
        Err(LaminateThermoelasticError::NotConverged { .. })
    ));
    assert_eq!(coupled.time_s(), 0.0);
    assert_eq!(
        coupled.state().shell().response().displacement_6n,
        initial_displacement
    );
    let mut checkpoints = 0;
    let step = coupled
        .advance_to(
            &thermal,
            0.5,
            &load,
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut workspace,
            || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    assert_eq!(coupled.time_s(), 0.5);
    assert_eq!(coupled.state().cell_potentials().len(), 3);
    let temperatures = coupled.state().shell().temperatures_k();
    assert!((temperatures[0] - temperatures[2]).abs() > 1.0);
    assert!(coupled.state().shell().response().displacement_6n[20].abs() > 1e-4);
    assert!(coupled.state().transverse_shear_energy_j() > 1e-12);
    assert!(step.report.transverse_shear_work_j.abs() > 1e-12);
    close(
        step.report.transverse_shear_work_j,
        step.report.transverse_shear_energy_change_j,
        1e-8,
    );
    close(
        step.thermal.energy_j() + field_internal(coupled.state())
            - thermal.energy_j()
            - initial_field,
        step.report.external_mechanical_work_j + step.report.external_heat_j,
        1e-7,
    );
    close(step.report.combined_energy_residual_j, 0.0, 1e-7);
    assert!(step.report.total_absolute_heat_defect_j <= controls().energy_tolerance_j);
    for cell in &step.report.cells {
        close(cell.free_energy_residual_j, 0.0, 1e-7);
    }
    assert_eq!(thermal.cells(), &original_cells);
    // The successful state cannot be paired with the stale original heat owner.
    let accepted = coupled.state().shell().response().displacement_6n.clone();
    assert!(matches!(
        coupled.advance_to(
            &thermal,
            1.0,
            &load,
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut workspace,
            || false
        ),
        Err(LaminateThermoelasticError::Invalid(_))
    ));
    assert_eq!(coupled.time_s(), 0.5);
    assert_eq!(coupled.state().shell().response().displacement_6n, accepted);
    // Repeat the same physical proposal from a fresh accepted owner and cancel
    // at its final checkpoint, after the thermal and mechanical candidates.
    let mut cancelled = ReciprocalLayeredThermalShell::try_new(
        &shell,
        &thermal,
        &[0.0; 24],
        &ranges,
        0.02,
        controls(),
        || false,
    )
    .unwrap();
    let mut retry = 0;
    assert!(matches!(
        cancelled.advance_to(
            &thermal,
            0.5,
            &load,
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut workspace,
            || {
                retry += 1;
                retry == checkpoints
            }
        ),
        Err(LaminateThermoelasticError::Cancelled)
    ));
    assert_eq!(cancelled.time_s(), 0.0);
    assert_eq!(
        cancelled.state().shell().response().displacement_6n,
        initial_displacement
    );
    assert_eq!(thermal.cells(), &original_cells);
}

#[test]
fn g0_reciprocal_admission_refuses_changing_elasticity_and_nonpositive_combined_capacity() {
    for excessive_expansion in [false, true] {
        let cards = [if excessive_expansion {
            card_with_coefficients("capacity refusal", 0.0, 1e-3, 0.0, CP)
        } else {
            card("variable elastic source")
        }];
        let body = prepare_body(&cards, &[1], 1.0);
        let laws = body.thermal_laws(false).unwrap();
        let bound = body.bind(&laws, ids()).unwrap();
        let thermal = initial(&bound, &[T0]);
        let shell = LayeredThermalShell::try_new(
            &bound,
            mechanical(&cards),
            ReferenceSurface::Midsurface,
            discretization(),
            limits(),
            || false,
        )
        .unwrap();
        let result = ReciprocalLayeredThermalShell::try_new(
            &shell,
            &thermal,
            &[0.0; 24],
            &[[280.0, 360.0]],
            0.1,
            controls(),
            || false,
        );
        let expected = if excessive_expansion {
            "capacity is not positive"
        } else {
            "fixed elastic sources"
        };
        assert!(
            matches!(result,Err(LaminateThermoelasticError::Invalid(message)) if message.contains(expected))
        );
        assert_eq!(thermal.time_s(), 0.0);
        close(
            bound.cell_temperatures(&thermal, || Ok(())).unwrap()[0],
            T0,
            0.0,
        );
    }
}
