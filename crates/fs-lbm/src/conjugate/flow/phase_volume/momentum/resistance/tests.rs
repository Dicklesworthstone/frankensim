//! G1/G3/G4: nonlinear physical drag heat, phase volume and one accepted impulse.

use super::super::tests::{close, controls as momentum_controls};
use super::{ResistiveThermalVolumeMomentumControls, ResistiveThermalVolumeMomentumStep};
use crate::conjugate::ChtError;
use crate::conjugate::domain::VoxelDomain;
use crate::conjugate::flow::{
    CoexistenceVolumeInput, FlowFace, FlowField, StaggeredBoundaryPressure,
    ThermalVolumeMomentumInput,
};
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};

fn curve() -> EquilibriumEnthalpyPhaseCurve {
    // The isothermal interval has v(h)=1/4+h/400, not affine density.
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([211; 32]),
        [
            (-100.0, 400.0, 0.0, 4.0),
            (0.0, 500.0, 0.0, 4.0),
            (100.0, 500.0, 1.0, 2.0),
            (200.0, 600.0, 1.0, 2.0),
        ]
        .into_iter()
        .map(|(h, temperature, liquid, density)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: temperature,
            liquid_mass_fraction: liquid,
            bulk_density_kg_m3: density,
        })
        .collect(),
    )
    .unwrap()
}

fn controls() -> ResistiveThermalVolumeMomentumControls {
    let mut thermal = momentum_controls();
    thermal.mechanics.projection.duration_s = 0.5;
    thermal.volume.thermal.volume.transport.projection = thermal.mechanics.projection;
    ResistiveThermalVolumeMomentumControls {
        thermal,
        max_coupling_iterations: 32,
        heat_relaxation: 1.0,
        maximum_cell_heat_defect_j: 1e-11,
    }
}

fn input<'a>(
    mass: &'a [f64],
    phases: &'a [Option<EquilibriumPhaseState>],
    increments: &'a [f64],
    momentum: &'a [[f64; 3]],
    force: &'a [[f64; 3]],
) -> ThermalVolumeMomentumInput<'a> {
    ThermalVolumeMomentumInput {
        volume: CoexistenceVolumeInput {
            mass_kg: mass,
            phase_states: phases,
            net_enthalpy_increment_j: increments,
            boundary_inflow: [None; 6],
        },
        momentum_kg_m_s: momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: force,
        boundary_velocity_m_s: [None; 6],
    }
}

fn positive_opening(axis: usize) -> [FlowFace; 6] {
    let mut faces = [FlowFace::Wall; 6];
    faces[2 * axis + 1] = FlowFace::Free;
    faces
}

fn physical_heat(external: f64) -> f64 {
    // One open half-dual has B=beta*V/2=10 kg/s. With dt=1/2 and
    // S=(E+H)/(400*dt), H=dt B S^2=(E+H)^2/8000. This is the stable
    // rationalized LOWER quadratic root, independent of the iteration.
    let c = 1.0 / 8000.0;
    2.0 * c * external * external / (1.0 - 2.0 * c * external + (1.0 - 4.0 * c * external).sqrt())
}

#[test]
fn drag_heat_and_phase_expansion_close_the_same_flow_and_inventory_on_every_axis() {
    let chart = curve();
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [16.0 / 5.0];
    let phases = [Some(chart.state_at_specific_enthalpy(25.0).unwrap())];
    let increments = [20.0];
    let zero = [[0.0; 3]];
    let heat = physical_heat(increments[0]);
    let flux = (increments[0] + heat) / 200.0;
    let final_mass = mass[0] / (1.0 + 0.5 * flux);
    let final_h = 25.0 + (increments[0] + heat) / mass[0];
    for axis in 0..3 {
        let result = field
            .advance_piecewise_thermal_volume_momentum_with_resistance(
                &domain,
                positive_opening(axis),
                &chart,
                input(&mass, &phases, &increments, &zero, &zero),
                &[20.0],
                controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        let volume = &result.thermal.volume;
        let transport = &volume.transport;
        let native = &result.resistance;
        assert!(result.coupling.coupling_iterations > 1);
        assert!(result.coupling.coupling_iterations <= controls().max_coupling_iterations);
        assert!(
            result.coupling.maximum_cell_heat_defect_j <= controls().maximum_cell_heat_defect_j
        );
        close(result.coupling.absolute_heat_coupling_defect_j, 0.0);
        assert_eq!(
            result.coupling.absolute_constitutive_impulse_defect_kg_m_s,
            0.0
        );
        assert_eq!(result.coupling.absolute_constitutive_heat_defect_j, 0.0);
        assert_eq!(result.resistance_pa_s_m2, [20.0]);
        close(result.caloric_resistance_heat_j[0], heat);
        close(native.cell_resistance_heat_j[0], heat);
        close(native.resistance_heat_j, heat);
        close(volume.report.requested_enthalpy_increment_j, 20.0 + heat);
        close(volume.report.represented_enthalpy_increment_j, 20.0 + heat);
        close(volume.volume_source_m3_s[0], flux);
        close(
            transport.projection.field.outward(&domain, 0, 2 * axis + 1),
            flux,
        );
        close(
            transport.projection.field.outward(&domain, 0, 2 * axis),
            0.0,
        );
        close(transport.mass_kg[0], final_mass);
        close(transport.specific_enthalpy_j_kg[0], final_h);
        close(transport.report.incoming_mass_kg, 0.0);
        close(transport.report.outgoing_mass_kg, mass[0] - final_mass);
        let endpoint = volume.phase_states[0].unwrap();
        close(endpoint.specific_enthalpy_j_kg(), final_h);
        close(endpoint.temperature_k(), 500.0);
        close(endpoint.liquid_mass_fraction(), final_h / 100.0);
        close(
            endpoint.bulk_density_kg_m3(),
            (0.25 + final_h / 400.0).recip(),
        );
        close(transport.mass_kg[0] / endpoint.bulk_density_kg_m3(), 1.0);
        assert_eq!(endpoint.phase_curve_identity(), chart.identity());
        close(
            volume.report.material_enthalpy_change_j,
            final_mass * final_h - mass[0] * 25.0,
        );
        close(
            volume.report.material_enthalpy_change_j,
            20.0 + heat - transport.report.outgoing_mass_kg * final_h,
        );

        // Original half-mass M/2=1.6 supplies inertia; dt B=5 supplies
        // resistance. Pressure contributes (1.6+5)Q and drag removes 5Q.
        // Only the accepted remaining impulse 1.6Q reaches original cell P.
        for component in 0..3 {
            let speed = if component == axis { flux } else { 0.0 };
            close(native.pressure_impulse_kg_m_s[component], 6.6 * speed);
            close(native.resistance_impulse_kg_m_s[component], -5.0 * speed);
            close(native.momentum_change_kg_m_s[component], 1.6 * speed);
            close(
                result.feedback.cell_momentum_change_kg_m_s[component],
                1.6 * speed,
            );
            close(
                result.momentum.momentum_kg_m_s[0][component],
                final_mass * speed / 2.0,
            );
        }
        close(native.kinetic_energy_change_j, 0.8 * flux * flux);
        close(native.backward_euler_loss_j, 0.8 * flux * flux);
        close(result.feedback.numerical_scatter_loss_j, 0.4 * flux * flux);
        close(
            volume.report.required_source_pressure_work_j,
            6.6 * flux * flux,
        );
        close(
            native.kinetic_energy_change_j + native.backward_euler_loss_j + heat,
            volume.report.required_source_pressure_work_j,
        );
        close(result.feedback.composed_energy_defect_j, 0.0);
        close(result.momentum.report.numerical_transport_loss_j, 0.0);
        assert_eq!(result.momentum.report.represented_impulse_kg_m_s, [0.0; 3]);
    }
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [16.0 / 5.0]);
    assert_eq!(
        phases,
        [Some(chart.state_at_specific_enthalpy(25.0).unwrap())]
    );
    assert_eq!(increments, [20.0]);
}

#[test]
fn zero_drag_reduces_to_original_thermal_momentum_with_large_reservoir_gauge() {
    let chart = curve();
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let faces = positive_opening(0);
    let mass = [16.0 / 5.0];
    let phases = [Some(chart.state_at_specific_enthalpy(25.0).unwrap())];
    let increments = [-20.0];
    let zero = [[0.0; 3]];
    let mut data = input(&mass, &phases, &increments, &zero, &zero);
    data.volume.boundary_inflow[1] = phases[0];
    data.boundary_velocity_m_s[1] = Some([-0.1, 0.0, 0.0]);
    let mut baseline: Option<ResistiveThermalVolumeMomentumStep> = None;
    for reference in [0.0, 1e12] {
        let x = [None, Some(reference)];
        let transverse = [None; 2];
        let profile = StaggeredBoundaryPressure {
            x_pa: &x,
            y_pa: &transverse,
            z_pa: &transverse,
        };
        let result = field
            .advance_piecewise_thermal_volume_momentum_with_resistance_and_pressure_profile(
                &domain,
                faces,
                profile,
                &chart,
                data,
                &[0.0],
                controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        let original = field
            .advance_piecewise_thermal_volume_momentum_with_pressure_profile(
                &domain,
                faces,
                profile,
                &chart,
                data,
                controls().thermal,
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        assert_eq!(result.thermal, original.thermal);
        assert_eq!(result.momentum, original.momentum);
        assert_eq!(result.reconstruction, original.reconstruction);
        assert_eq!(result.coupling.coupling_iterations, 1);
        assert_eq!(result.caloric_resistance_heat_j, [0.0]);
        assert_eq!(result.resistance.resistance_heat_j, 0.0);
        assert_eq!(result.coupling.absolute_heat_coupling_defect_j, 0.0);
        let volume = &result.thermal.volume;
        close(volume.volume_source_m3_s[0], -0.1);
        close(volume.transport.mass_kg[0], 84.0 / 25.0);
        close(volume.transport.specific_enthalpy_j_kg[0], 400.0 / 21.0);
        close(volume.transport.report.incoming_mass_kg, 4.0 / 25.0);
        close(result.momentum.momentum_kg_m_s[0][0], -22.0 / 125.0);
        if let Some(prior) = &baseline {
            assert_eq!(
                volume.transport.projection.field,
                prior.thermal.volume.transport.projection.field
            );
            assert_eq!(volume.phase_states, prior.thermal.volume.phase_states);
            assert_eq!(result.momentum, prior.momentum);
            assert_eq!(result.resistance, prior.resistance);
            assert_eq!(result.feedback, prior.feedback);
        } else {
            baseline = Some(result);
        }
    }
}

#[test]
fn bounded_heat_solve_and_late_incoming_momentum_failure_preserve_retry() {
    let chart = curve();
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let original_field = field.clone();
    let faces = positive_opening(0);
    let mass = [16.0 / 5.0];
    let phases = [Some(chart.state_at_specific_enthalpy(25.0).unwrap())];
    let increments = [-20.0];
    let zero = [[0.0; 3]];
    let mut data = input(&mass, &phases, &increments, &zero, &zero);
    data.volume.boundary_inflow[1] = phases[0];
    data.boundary_velocity_m_s[1] = Some([-0.1, 0.0, 0.0]);
    let run = |candidate: ThermalVolumeMomentumInput<'_>,
               control: ResistiveThermalVolumeMomentumControls,
               gate: &CancelGate| {
        field.advance_piecewise_thermal_volume_momentum_with_resistance(
            &domain,
            faces,
            &chart,
            candidate,
            &[20.0],
            control,
            gate,
        )
    };
    let expected = run(data, controls(), &CancelGate::new_clock_free()).unwrap();
    close(expected.resistance.resistance_heat_j, physical_heat(-20.0));
    let mut exhausted = controls();
    exhausted.max_coupling_iterations = 1;
    assert!(run(data, exhausted, &CancelGate::new_clock_free()).is_err());
    let mut invalid = controls();
    invalid.heat_relaxation = 0.0;
    assert!(run(data, invalid, &CancelGate::new_clock_free()).is_err());
    let mut mismatched = controls();
    mismatched.thermal.mechanics.projection.duration_s *= 2.0;
    assert!(run(data, mismatched, &CancelGate::new_clock_free()).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        run(data, controls(), &cancelled),
        Err(ChtError::Cancelled)
    ));

    // Thermal/volume and native feedback can accept this original inventory;
    // only the actual incoming mass route reveals the missing velocity.
    let mut missing_velocity = data;
    missing_velocity.boundary_velocity_m_s = [None; 6];
    assert!(matches!(
        run(missing_velocity, controls(), &CancelGate::new_clock_free()),
        Err(ChtError::InvalidInput { field: "flow.material_advection", reason })
            if reason.contains("incoming mass route lacks an authored velocity")
    ));
    assert_eq!(field, original_field);
    assert_eq!(mass, [16.0 / 5.0]);
    assert_eq!(
        phases,
        [Some(chart.state_at_specific_enthalpy(25.0).unwrap())]
    );
    assert_eq!(increments, [-20.0]);
    assert_eq!(zero, [[0.0; 3]]);
    assert_eq!(
        run(data, controls(), &CancelGate::new_clock_free()).unwrap(),
        expected,
    );
}
