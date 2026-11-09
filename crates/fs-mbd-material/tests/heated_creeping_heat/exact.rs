//! Same-source conduction, exact free-spin relaxation and atomic observers.
use super::*;
use fs_mbd_material::thermal::heating::motion::creeping_flow::HeatedExactCreepingControls;

fn exact_controls() -> HeatedExactCreepingControls {
    HeatedExactCreepingControls {
        maximum_step_energy_defect_j: 1e-11,
    }
}

#[test]
fn stiff_spin_and_original_conduction_share_the_finite_cell_without_double_heat() {
    let material = card(CP, 0.0, 1000.0);
    let chart = caloric(&material);
    let solid = transport(&material, 1e5, 360.0);
    let reference = reference(&material, 300.0, RADIUS);
    let medium = ThermalFluidTransportState::try_new(
        FluidTransportState::try_new(0.001, 1.0, 10.0, 1000.0).unwrap(),
        320.0,
        None,
    )
    .unwrap();
    let sphere =
        CreepingSphere::new(&medium, RADIUS, Gravity::ZERO, Vec3::ZERO, flow_limits()).unwrap();
    let dt = 0.02;
    let original = body(
        &reference,
        &material,
        &chart,
        CP,
        300.0,
        RADIUS,
        dt,
        Gravity::ZERO,
        0.2,
    );
    assert!(
        original
            .advance_creeping_sphere_in_thermal_fluid(
                &sphere,
                exchange(&solid),
                flow_controls(),
                dt,
                config(dt),
                || Ok(())
            )
            .is_err()
    );
    let (next, report) = original
        .advance_creeping_sphere_in_thermal_fluid_exact(
            &sphere,
            exchange(&solid),
            exact_controls(),
            dt,
            config(dt),
            || Ok(()),
        )
        .unwrap();
    let inertia = original
        .principal_frame()
        .properties()
        .principal_inertia_body()
        .z;
    let resistance = 8.0 * std::f64::consts::PI * 1000.0 * RADIUS.powi(3);
    let expected_spin = 0.2 * (-dt * resistance / inertia).exp();
    let actual_spin = next
        .principal_frame()
        .properties()
        .angular_velocity_body_checked(next.state().angular_momentum_body())
        .unwrap()
        .z;
    close(actual_spin, expected_spin, 1e-14);
    let conductance = 4.0 * std::f64::consts::PI * 10.0 * RADIUS;
    let capacity = reference.inertia().mass() * CP;
    // The original conduction observer adds its inverse stable time to the
    // mechanical cap. Constant Cp, k and radius make the two actual Euler
    // intervals known analytically; mechanical relaxation composes exactly.
    let first = 1.0 / (1.0 / dt + conductance / capacity / 0.8);
    let last = dt - first;
    assert_eq!(report.coupled.thermal.substeps, 2);
    let heat = conductance * 20.0 * first
        + conductance * 20.0 * (1.0 - conductance * first / capacity) * last;
    close(report.reservoir_heat_into_body_j, heat, 1e-10);
    close(report.coupled.thermal.external_heat_j, heat, 1e-10);
    close(
        next.thermal().temperature_k(CellId(0)).unwrap(),
        300.0 + heat / (reference.inertia().mass() * CP),
        1e-11,
    );
    close(
        report.coupled.viscous_heat_j,
        0.5 * inertia * 0.2_f64.powi(2) * -(-2.0 * dt * resistance / inertia).exp_m1(),
        1e-12,
    );
    assert!(report.coupled.viscous_heat_j > 0.0);
    close(report.coupled.combined_residual_j, 0.0, 1e-8);
    assert_eq!(
        next.configuration().inertia().mass(),
        reference.inertia().mass()
    );
}

#[test]
fn exact_conduction_composition_retains_expansion_sources_and_atomic_failure() {
    let material = card(CP, 0.002, 1000.0);
    let chart = caloric(&material);
    let solid = transport(&material, 1e5, 360.0);
    let reference = reference(&material, 300.0, RADIUS);
    let medium = ThermalFluidTransportState::try_new(
        FluidTransportState::try_new(0.001, 1.0, 10.0, 1000.0).unwrap(),
        320.0,
        None,
    )
    .unwrap();
    let sphere =
        CreepingSphere::new(&medium, RADIUS, Gravity::ZERO, Vec3::ZERO, flow_limits()).unwrap();
    let dt = 0.02;
    let original = body(
        &reference,
        &material,
        &chart,
        CP,
        300.0,
        RADIUS,
        dt,
        Gravity::ZERO,
        0.2,
    );
    let heater = Heater {
        watts: 2.0,
        steps: 0,
        refuse_at: usize::MAX,
    };
    let (next, generated, report) = original
        .advance_creeping_sphere_in_thermal_fluid_exact_with(
            &heater,
            &sphere,
            exchange(&solid),
            exact_controls(),
            3.0 * dt,
            config(dt),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(generated.steps, report.coupled.thermal.substeps);
    assert!(generated.steps > 3);
    close(
        report.coupled.thermal.external_heat_j,
        report.reservoir_heat_into_body_j + 2.0 * 3.0 * dt,
        1e-9,
    );
    let temperature = next.thermal().temperature_k(CellId(0)).unwrap();
    close(
        report.after.radius_m,
        RADIUS * (1.0 + 0.002 * (temperature - 300.0)),
        1e-13,
    );
    close(
        report.after.conductance_w_k / report.before.conductance_w_k,
        report.after.radius_m / RADIUS,
        1e-13,
    );
    close(report.coupled.combined_residual_j, 0.0, 1e-8);
    let failure = Heater {
        refuse_at: 2,
        ..heater.clone()
    };
    assert!(
        original
            .advance_creeping_sphere_in_thermal_fluid_exact_with(
                &failure,
                &sphere,
                exchange(&solid),
                exact_controls(),
                3.0 * dt,
                config(dt),
                || Ok(())
            )
            .is_err()
    );
    assert_eq!(failure.steps, 0);
    assert_eq!(original.thermal().time_s(), 0.0);
    assert_eq!(original.thermal().temperature_k(CellId(0)).unwrap(), 300.0);
}
