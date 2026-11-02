//! G1/G3/G4: actual phase-volume flow, conservative inventories and refusal.
use super::super::DensityProjectionControls;
use super::*;
use fs_blake3::ContentHash;
use fs_material::phase::EnthalpyPhaseKnot;

fn chart(offset: f64, liquid_density: f64, identity: u8) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([identity; 32]),
        [
            (-2.0, 300.0, 0.0, 2.0),
            (0.0, 400.0, 0.0, 2.0),
            (4.0, 400.0, 1.0, liquid_density),
            (6.0, 500.0, 1.0, liquid_density),
        ]
        .into_iter()
        .map(|(h, t, liquid, rho)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: offset + h,
            temperature_k: t,
            liquid_mass_fraction: liquid,
            bulk_density_kg_m3: rho,
        })
        .collect(),
    )
    .unwrap()
}

fn controls() -> CoexistenceVolumeControls {
    CoexistenceVolumeControls {
        transport: MaterialAdvectionControls {
            projection: DensityProjectionControls {
                duration_s: 0.5,
                relative_tolerance: 1e-12,
                maximum_volume_defect_m3_s: 1e-11,
                maximum_energy_defect_j: 1e-10,
            },
            max_cells: 100,
            max_connections: 1000,
            max_iterations: 1000,
            relative_tolerance: 1e-12,
            maximum_mass_defect_kg: 1e-10,
            maximum_enthalpy_defect_j: 1e-9,
        },
        maximum_cell_volume_defect_m3: 1e-10,
    }
}

fn input<'a>(
    mass: &'a [f64],
    states: &'a [Option<EquilibriumPhaseState>],
    increments: &'a [f64],
) -> CoexistenceVolumeInput<'a> {
    CoexistenceVolumeInput {
        mass_kg: mass,
        phase_states: states,
        net_enthalpy_increment_j: increments,
        boundary_inflow: [None; 6],
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-11 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn outlet(axis: usize) -> [FlowFace; 6] {
    let mut faces = [FlowFace::Wall; 6];
    faces[2 * axis + 1] = FlowFace::Free;
    faces
}

#[test]
fn melting_exports_original_material_and_continuation_refills_the_same_phase_inventory() {
    // v(h)=1/2+h/8. M0=4/3 fills V=1 at h0=2. The prescribed
    // delta_H=M0 gives h*=3, delta_V=1/6. The original donor equation
    // yields M1=(4/3)/(1+1/6)=8/7, h1=3, and exported mass 4/21.
    let curve = chart(0.0, 1.0, 71);
    for axis in 0..3 {
        let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
        let field = FlowField::quiescent(&domain);
        let original = curve.state_at_specific_enthalpy(2.0).unwrap();
        let mass = [4.0 / 3.0];
        let states = [Some(original)];
        let increments = [mass[0]];
        let hot = field
            .advance_coexistence_volume(
                &domain,
                outlet(axis),
                &curve,
                1,
                input(&mass, &states, &increments),
                controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        close(hot.volume_source_m3_s[0], 1.0 / 3.0);
        close(
            hot.transport
                .projection
                .field
                .outward(&domain, 0, 2 * axis + 1),
            1.0 / 3.0,
        );
        close(hot.transport.mass_kg[0], 8.0 / 7.0);
        close(hot.transport.specific_enthalpy_j_kg[0], 3.0);
        close(hot.transport.report.outgoing_mass_kg, 4.0 / 21.0);
        close(hot.transport.report.incoming_mass_kg, 0.0);
        close(hot.transport.report.carried_enthalpy_j, -4.0 / 7.0);
        close(hot.report.material_enthalpy_change_j, 16.0 / 21.0);
        close(hot.report.requested_enthalpy_increment_j, 4.0 / 3.0);
        close(hot.report.represented_enthalpy_increment_j, 4.0 / 3.0);
        close(hot.report.required_source_pressure_work_j, 2.0 / 27.0);
        assert!(hot.transport.projection.report.projection_loss_j > 0.0);
        let phase = hot.phase_states[0].unwrap();
        close(phase.temperature_k(), 400.0);
        close(phase.liquid_mass_fraction(), 0.75);
        assert_eq!(phase.phase_curve_identity(), curve.identity());
        close(hot.transport.mass_kg[0] / phase.bulk_density_kg_m3(), 1.0);
        close(hot.report.absolute_enthalpy_defect_j, 0.0);
        assert_eq!(field, FlowField::quiescent(&domain));
        assert_eq!(states, [Some(original)]);

        // Continue from the returned masses and phases, cooling to h*=2.
        // The same free face becomes an inlet with the original h=2 state.
        let cooling = [-hot.transport.mass_kg[0]];
        let mut next_input = input(&hot.transport.mass_kg, &hot.phase_states, &cooling);
        next_input.boundary_inflow[2 * axis + 1] = Some(original);
        let cold = hot
            .transport
            .projection
            .field
            .advance_coexistence_volume(
                &domain,
                outlet(axis),
                &curve,
                1,
                next_input,
                controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        close(cold.volume_source_m3_s[0], -2.0 / 7.0);
        close(cold.transport.mass_kg[0], mass[0]);
        close(cold.transport.specific_enthalpy_j_kg[0], 2.0);
        close(cold.transport.report.incoming_mass_kg, 4.0 / 21.0);
        close(cold.transport.report.outgoing_mass_kg, 0.0);
        close(cold.phase_states[0].unwrap().liquid_mass_fraction(), 0.5);
        close(cold.report.maximum_final_volume_defect_m3, 0.0);
    }
}

#[test]
fn contracting_mixture_admits_same_chart_inflow_and_mixes_extensive_enthalpy() {
    let curve = chart(0.0, 1.0, 72);
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let mass = [4.0 / 3.0];
    let states = [Some(curve.state_at_specific_enthalpy(2.0).unwrap())];
    let delta = [-4.0 / 3.0];
    let incoming = curve.state_at_specific_enthalpy(3.0).unwrap();
    let mut data = input(&mass, &states, &delta);
    data.boundary_inflow[1] = Some(incoming);
    let result = FlowField::quiescent(&domain)
        .advance_coexistence_volume(
            &domain,
            outlet(0),
            &curve,
            1,
            data,
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // h*=1 gives delta_V=-1/6. Inflow rho(3)=8/7 adds 4/21 kg.
    // M1=32/21, H1=4/3+3*(4/21)=40/21, hence h1=5/4.
    close(result.volume_source_m3_s[0], -1.0 / 3.0);
    close(result.transport.mass_kg[0], 32.0 / 21.0);
    close(result.transport.specific_enthalpy_j_kg[0], 1.25);
    close(result.transport.report.incoming_mass_kg, 4.0 / 21.0);
    close(result.transport.report.carried_enthalpy_j, 4.0 / 7.0);
    close(result.report.material_enthalpy_change_j, -16.0 / 21.0);
    close(
        result.phase_states[0].unwrap().bulk_density_kg_m3(),
        32.0 / 21.0,
    );
    close(result.report.maximum_final_volume_defect_m3, 0.0);
}

#[test]
fn balanced_sealed_melting_and_freezing_conserve_total_mass_and_enthalpy() {
    let curve = chart(0.0, 1.0, 73);
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let mass = [4.0 / 3.0; 2];
    let states = [Some(curve.state_at_specific_enthalpy(2.0).unwrap()); 2];
    let delta = [4.0 / 3.0, -4.0 / 3.0];
    let result = FlowField::quiescent(&domain)
        .advance_coexistence_volume(
            &domain,
            [FlowFace::Wall; 6],
            &curve,
            1,
            input(&mass, &states, &delta),
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    close(result.transport.projection.field.flux_x(1, 0, 0), 1.0 / 3.0);
    close(result.transport.mass_kg[0], 8.0 / 7.0);
    close(result.transport.mass_kg[1], 32.0 / 21.0);
    close(result.transport.specific_enthalpy_j_kg[0], 3.0);
    close(result.transport.specific_enthalpy_j_kg[1], 1.25);
    close(result.transport.mass_kg.iter().sum(), 8.0 / 3.0);
    close(
        result
            .transport
            .mass_kg
            .iter()
            .zip(&result.transport.specific_enthalpy_j_kg)
            .map(|(m, h)| m * h)
            .sum(),
        16.0 / 3.0,
    );
    close(result.transport.report.incoming_mass_kg, 0.0);
    close(result.transport.report.outgoing_mass_kg, 0.0);
    close(result.report.material_enthalpy_change_j, 0.0);
    close(result.report.required_source_pressure_work_j, 4.0 / 27.0);
    for (i, phase) in result.phase_states.iter().enumerate() {
        let phase = phase.unwrap();
        close(
            result.transport.mass_kg[i] / phase.bulk_density_kg_m3(),
            1.0,
        );
        close(phase.temperature_k(), 400.0);
        assert_eq!(phase.phase_curve_identity(), curve.identity());
    }
}

#[test]
fn equal_phase_densities_preserve_the_original_material_transport_exactly() {
    let curve = chart(0.0, 2.0, 74);
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let faces = [
        FlowFace::Fixed,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let field =
        FlowField::from_face_velocity(&domain, faces, |axis, _| if axis == 0 { 0.25 } else { 0.0 })
            .unwrap();
    let mass = [2.0];
    let states = [Some(curve.state_at_specific_enthalpy(2.0).unwrap())];
    let delta = [2.0];
    let incoming = curve.state_at_specific_enthalpy(3.0).unwrap();
    let mut data = input(&mass, &states, &delta);
    data.boundary_inflow[0] = Some(incoming);
    let phase = field
        .advance_coexistence_volume(
            &domain,
            faces,
            &curve,
            1,
            data,
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let mut inlets = [None; 6];
    inlets[0] = Some(AdvectedMaterial {
        density_kg_m3: 2.0,
        specific_enthalpy_j_kg: 3.0,
    });
    let native = field
        .project_and_advect_material(
            &domain,
            faces,
            MaterialAdvectionInput {
                mass_kg: &mass,
                specific_enthalpy_j_kg: &[3.0],
                volume_sources: &[MaterialVolumeSource::None],
                boundary_inflow: inlets,
            },
            controls().transport,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(phase.transport, native);
    assert_eq!(phase.volume_source_m3_s, [0.0]);
    close(phase.phase_states[0].unwrap().liquid_mass_fraction(), 0.75);
}

#[test]
fn admitted_volume_roundoff_does_not_drive_a_sealed_zero_input_continuation() {
    let curve = chart(0.0, 1.0, 79);
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let original = curve.state_at_specific_enthalpy(2.0).unwrap();
    // An admitted endpoint may carry one representable mass increment of EOS
    // roundoff. No caloric change means no physical expansion or contraction.
    let mass = [original.bulk_density_kg_m3().next_up()];
    let states = [Some(original)];
    let field = FlowField::quiescent(&domain);
    let first = field
        .advance_coexistence_volume(
            &domain,
            [FlowFace::Wall; 6],
            &curve,
            1,
            input(&mass, &states, &[0.0]),
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(first.volume_source_m3_s, [0.0]);
    assert_eq!(first.transport.mass_kg, mass);
    assert_eq!(first.phase_states, states);
    assert!(first.report.maximum_initial_volume_defect_m3 > 0.0);
    assert_eq!(
        first.report.maximum_initial_volume_defect_m3,
        first.report.maximum_final_volume_defect_m3
    );
    assert_eq!(first.report.required_source_pressure_work_j, 0.0);
    let continued = first
        .transport
        .projection
        .field
        .advance_coexistence_volume(
            &domain,
            [FlowFace::Wall; 6],
            &curve,
            1,
            input(&first.transport.mass_kg, &first.phase_states, &[0.0]),
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(continued, first);
}

#[test]
fn caloric_datum_shift_preserves_phase_volume_mass_and_pressure() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let results: Vec<_> = [0.0, (1_u64 << 32) as f64]
        .into_iter()
        .map(|offset| {
            let curve = chart(offset, 1.0, 75);
            let mut controls = controls();
            controls.transport.maximum_enthalpy_defect_j = 1e-5;
            FlowField::quiescent(&domain)
                .advance_coexistence_volume(
                    &domain,
                    outlet(0),
                    &curve,
                    1,
                    input(
                        &[4.0 / 3.0],
                        &[Some(
                            curve.state_at_specific_enthalpy(offset + 2.0).unwrap(),
                        )],
                        &[4.0 / 3.0],
                    ),
                    controls,
                    &CancelGate::new_clock_free(),
                )
                .unwrap()
        })
        .collect();
    let a = &results[0];
    let b = &results[1];
    assert_eq!(a.transport.mass_kg, b.transport.mass_kg);
    assert_eq!(a.transport.projection, b.transport.projection);
    assert_eq!(a.volume_source_m3_s, b.volume_source_m3_s);
    assert_eq!(
        a.phase_states[0].unwrap().liquid_mass_fraction(),
        b.phase_states[0].unwrap().liquid_mass_fraction()
    );
    assert_eq!(
        a.report.required_source_pressure_work_j,
        b.report.required_source_pressure_work_j
    );
    close(
        b.transport.specific_enthalpy_j_kg[0] - (1_u64 << 32) as f64,
        3.0,
    );
}

#[test]
fn phase_domain_identity_inventory_caloric_and_cancellation_refusals_are_atomic() {
    let curve = chart(0.0, 1.0, 76);
    let other = chart(0.0, 1.0, 77);
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [4.0 / 3.0];
    let states = [Some(curve.state_at_specific_enthalpy(2.0).unwrap())];
    let delta = [4.0 / 3.0];
    let data = input(&mass, &states, &delta);
    let run = |interval, faces, data, controls, gate: &CancelGate| {
        field.advance_coexistence_volume(&domain, faces, &curve, interval, data, controls, gate)
    };
    let gate = CancelGate::new_clock_free();
    let baseline = run(1, outlet(0), data, controls(), &gate).unwrap();
    for interval in [0, 2, 3, usize::MAX] {
        assert!(run(interval, outlet(0), data, controls(), &gate).is_err());
    }
    assert!(run(1, [FlowFace::Wall; 6], data, controls(), &gate).is_err());
    assert!(
        run(
            1,
            outlet(0),
            input(&[1.0], &states, &delta),
            controls(),
            &gate
        )
        .is_err()
    );
    assert!(
        run(
            1,
            outlet(0),
            input(&mass, &states, &[8.0]),
            controls(),
            &gate
        )
        .is_err()
    );
    let lost_increment = [f64::from_bits(1)];
    assert!(
        run(
            1,
            outlet(0),
            input(&mass, &states, &lost_increment),
            controls(),
            &gate
        )
        .is_err()
    );
    let alien = [Some(other.state_at_specific_enthalpy(2.0).unwrap())];
    assert!(
        run(
            1,
            outlet(0),
            input(&mass, &alien, &delta),
            controls(),
            &gate
        )
        .is_err()
    );
    let mut boundary = data;
    boundary.boundary_inflow[1] = alien[0];
    assert!(run(1, outlet(0), boundary, controls(), &gate).is_err());
    let mut strict = controls();
    strict.transport.maximum_enthalpy_defect_j = 0.0;
    let rounded_increment = [mass[0] / 10.0];
    assert!(
        run(
            1,
            outlet(0),
            input(&mass, &states, &rounded_increment),
            strict,
            &gate
        )
        .is_err()
    );
    assert_eq!(
        run(1, outlet(0), data, controls(), &gate).unwrap(),
        baseline
    );
    gate.request();
    assert!(run(1, outlet(0), data, controls(), &gate).is_err());
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [4.0 / 3.0]);
    assert_eq!(states[0].unwrap().specific_enthalpy_j_kg(), 2.0);
}

#[test]
fn independent_final_volume_gate_rejects_an_inaccurate_material_solve() {
    let curve = chart(0.0, 1.0, 78);
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let mass = vec![4.0 / 3.0; 9];
    let states = vec![Some(curve.state_at_specific_enthalpy(2.0).unwrap()); 9];
    let delta: Vec<_> = (0..9).map(|i| mass[i] * (0.1 + 0.1 * i as f64)).collect();
    let mut loose = controls();
    loose.transport.relative_tolerance = 0.5;
    loose.transport.maximum_mass_defect_kg = 100.0;
    loose.transport.maximum_enthalpy_defect_j = 100.0;
    let error = FlowField::quiescent(&domain)
        .advance_coexistence_volume(
            &domain,
            [FlowFace::Free; 6],
            &curve,
            1,
            input(&mass, &states, &delta),
            loose,
            &CancelGate::new_clock_free(),
        )
        .unwrap_err();
    assert!(
        matches!(error,ChtError::InvalidInput {field:"flow.coexistence_volume",ref reason}
        if reason.contains("transported mass and phase density")),
        "{error:?}"
    );
}

fn sensible_chart(offset: f64) -> EquilibriumEnthalpyPhaseCurve {
    use fs_material::phase::{HeatCapacityKnot, SolidLiquidPhase};
    EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
        ContentHash([91; 32]),
        SolidLiquidPhase::Liquid,
        offset,
        &[
            HeatCapacityKnot {
                temperature_k: 300.0,
                specific_heat_capacity_j_kg_k: 1.0,
                bulk_density_kg_m3: 2.0,
            },
            HeatCapacityKnot {
                temperature_k: 304.0,
                specific_heat_capacity_j_kg_k: 1.0,
                bulk_density_kg_m3: 1.0,
            },
        ],
        1e-12,
    )
    .unwrap()
}

fn thermal_controls() -> ThermalVolumeControls {
    ThermalVolumeControls {
        volume: controls(),
        max_volume_iterations: 100,
        source_relaxation: 1.0,
    }
}

#[test]
fn sensible_heat_and_cooling_drive_original_density_volume_on_every_axis() {
    let curve = sensible_chart(0.0);
    assert_eq!(
        curve.density_law_on_interval(0).unwrap(),
        PhaseIntervalDensityLaw::AffineDensity
    );
    assert!(curve.density_law_on_interval(1).is_err());
    assert!(curve.density_law_on_interval(usize::MAX).is_err());
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let mass = [2.0];
    let states = [Some(curve.state_at_specific_enthalpy(0.0).unwrap())];
    for axis in 0..3 {
        let field = FlowField::quiescent(&domain);
        let hot = field
            .advance_thermal_volume(
                &domain,
                outlet(axis),
                &curve,
                0,
                input(&mass, &states, &[2.0]),
                thermal_controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        // Cp=1 gives T=301, and the ORIGINAL rho(T) is 2-1/4.
        // The same M=2 occupies 8/7 before exporting 1/4 kg.
        let step = &hot.volume;
        close(step.transport.mass_kg[0], 7.0 / 4.0);
        close(step.transport.specific_enthalpy_j_kg[0], 1.0);
        close(step.phase_states[0].unwrap().temperature_k(), 301.0);
        close(
            step.phase_states[0].unwrap().bulk_density_kg_m3(),
            7.0 / 4.0,
        );
        close(step.volume_source_m3_s[0], 2.0 / 7.0);
        close(step.transport.report.outgoing_mass_kg, 1.0 / 4.0);
        close(step.report.material_enthalpy_change_j, 7.0 / 4.0);
        assert_eq!(hot.source_iterations, 1);
        assert_eq!(hot.mixing_volume_source_m3_s, [0.0]);
        assert_eq!(
            step.phase_states[0].unwrap().phase_curve_identity(),
            curve.identity()
        );
        let cooling = [-step.transport.mass_kg[0]];
        let mut data = input(&step.transport.mass_kg, &step.phase_states, &cooling);
        data.boundary_inflow[2 * axis + 1] = states[0];
        let cold = step
            .transport
            .projection
            .field
            .advance_thermal_volume(
                &domain,
                outlet(axis),
                &curve,
                0,
                data,
                thermal_controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        close(cold.volume.transport.mass_kg[0], mass[0]);
        close(cold.volume.transport.specific_enthalpy_j_kg[0], 0.0);
        close(cold.volume.transport.report.incoming_mass_kg, 1.0 / 4.0);
        close(cold.volume.volume_source_m3_s[0], -1.0 / 4.0);
        assert_eq!(cold.source_iterations, 1);
        assert_eq!(cold.mixing_volume_source_m3_s, [0.0]);
        assert_eq!(field, FlowField::quiescent(&domain));
    }
}

#[test]
fn sensible_inflow_mixing_solves_original_density_and_energy_quadratic() {
    let curve = sensible_chart(0.0);
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let mass = [1.5];
    let states = [Some(curve.state_at_specific_enthalpy(2.0).unwrap())];
    let increments = [-1.5];
    let mut data = input(&mass, &states, &increments);
    data.boundary_inflow[1] = Some(curve.state_at_specific_enthalpy(3.0).unwrap());
    let mut ctl = thermal_controls();
    ctl.volume.maximum_cell_volume_defect_m3 = 1e-12;
    let next = FlowField::quiescent(&domain)
        .advance_thermal_volume(
            &domain,
            outlet(0),
            &curve,
            0,
            data,
            ctl,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // With incoming volume x, M=1.5+1.25x and H=1.5+3.75x.
    // EOS M=2-H/(4M) gives M^2-1.25M-.75=0 independently
    // of the pressure/source iteration. Heat is supplied once, not per trial.
    let expected_mass = (5.0 + 73.0_f64.sqrt()) / 8.0;
    let incoming_volume = (expected_mass - 1.5) / 1.25;
    let expected_h = 3.0 - 3.0 / expected_mass;
    let step = &next.volume;
    close(step.transport.mass_kg[0], expected_mass);
    close(step.transport.specific_enthalpy_j_kg[0], expected_h);
    close(
        step.transport.report.incoming_mass_kg,
        expected_mass - mass[0],
    );
    close(step.volume_source_m3_s[0], -incoming_volume / 0.5);
    close(
        next.mixing_volume_source_m3_s[0],
        (-incoming_volume + 1.0 / 7.0) / 0.5,
    );
    assert!(next.source_iterations > 1 && next.source_iterations < 20);
    close(step.report.requested_enthalpy_increment_j, -1.5);
    close(step.report.represented_enthalpy_increment_j, -1.5);
    close(
        step.report.material_enthalpy_change_j,
        -1.5 + 3.0 * (expected_mass - mass[0]),
    );
    close(
        step.report.required_source_pressure_work_j,
        step.transport.projection.report.pressure_source_work_j,
    );
    assert!(step.report.maximum_final_volume_defect_m3 <= 1e-12);
    assert_eq!(mass, [1.5]);
    assert_eq!(states[0].unwrap().specific_enthalpy_j_kg(), 2.0);
}

#[test]
fn graded_sensible_mixing_preserves_internal_mass_and_enthalpy_on_every_axis() {
    let curve = sensible_chart(0.0);
    for axis in 0..3 {
        let mut widths = [vec![2.0], vec![3.0], vec![2.0]];
        widths[axis] = vec![1.0, 2.0];
        widths[(axis + 1) % 3] = vec![2.0];
        widths[(axis + 2) % 3] = vec![3.0];
        let domain = VoxelDomain::graded(widths).unwrap();
        let mut faces = outlet(axis);
        faces[2 * axis] = FlowFace::Fixed;
        let field =
            FlowField::from_face_velocity(&domain, faces, |a, _| if a == axis { 1.0 } else { 0.0 })
                .unwrap();
        let mass = [12.0, 18.0];
        let states = [
            Some(curve.state_at_specific_enthalpy(0.0).unwrap()),
            Some(curve.state_at_specific_enthalpy(2.0).unwrap()),
        ];
        let mut data = input(&mass, &states, &[0.0; 2]);
        data.boundary_inflow[2 * axis] = Some(curve.state_at_specific_enthalpy(3.0).unwrap());
        let mut ctl = thermal_controls();
        ctl.volume.maximum_cell_volume_defect_m3 = 1e-11;
        let next = field
            .advance_thermal_volume(
                &domain,
                faces,
                &curve,
                0,
                data,
                ctl,
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        // A=6, dt=.5: the fixed inlet supplies m_in=3*1.25.
        // Cell donor h is its old+incoming extensive enthalpy / total
        // available mass. Original EOS fixes retained M; the remainder exits.
        let incoming_mass = 3.75;
        let h0 = incoming_mass * 3.0 / (mass[0] + incoming_mass);
        let rho0 = 2.0 - h0 / 4.0;
        let m0 = 6.0 * rho0;
        let transferred = mass[0] + incoming_mass - m0;
        let h1 = (mass[1] * 2.0 + transferred * h0) / (mass[1] + transferred);
        let rho1 = 2.0 - h1 / 4.0;
        let m1 = 12.0 * rho1;
        let outgoing_mass = mass[1] + transferred - m1;
        let step = &next.volume;
        close(step.transport.mass_kg[0], m0);
        close(step.transport.mass_kg[1], m1);
        close(step.transport.specific_enthalpy_j_kg[0], h0);
        close(step.transport.specific_enthalpy_j_kg[1], h1);
        close(step.transport.report.incoming_mass_kg, incoming_mass);
        close(step.transport.report.outgoing_mass_kg, outgoing_mass);
        close(
            step.transport
                .projection
                .field
                .outward(&domain, 0, 2 * axis + 1)
                * 0.5
                * rho0,
            transferred,
        );
        close(
            step.report.material_enthalpy_change_j,
            incoming_mass * 3.0 - outgoing_mass * h1,
        );
        close(
            step.transport.report.mass_after_kg - step.transport.report.mass_before_kg,
            incoming_mass - outgoing_mass,
        );
        assert!(next.mixing_volume_source_m3_s.iter().all(|s| *s < 0.0));
        assert!(next.source_iterations > 1);
        assert_eq!(step.report.requested_enthalpy_increment_j, 0.0);
    }
}

#[test]
fn affine_sensible_and_coexistence_use_the_same_exact_volume_owner() {
    use fs_material::phase::SolidLiquidPhase;
    let liquid = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([92; 32]),
        SolidLiquidPhase::Liquid,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 300.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 2.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 4.0,
                temperature_k: 304.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1.0,
            },
        ],
    )
    .unwrap();
    assert_eq!(
        liquid.density_law_on_interval(0).unwrap(),
        PhaseIntervalDensityLaw::AffineSpecificVolume
    );
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let sensible = field
        .advance_thermal_volume(
            &domain,
            outlet(0),
            &liquid,
            0,
            input(
                &[2.0],
                &[Some(liquid.state_at_specific_enthalpy(0.0).unwrap())],
                &[2.0],
            ),
            thermal_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    close(sensible.volume.transport.mass_kg[0], 1.6);
    close(sensible.volume.volume_source_m3_s[0], 0.5);
    assert_eq!(sensible.source_iterations, 1);
    assert_eq!(sensible.mixing_volume_source_m3_s, [0.0]);

    let mixture = chart(0.0, 1.0, 93);
    let mass = [4.0 / 3.0];
    let states = [Some(mixture.state_at_specific_enthalpy(2.0).unwrap())];
    let increment = [mass[0]];
    let data = input(&mass, &states, &increment);
    let original = field
        .advance_coexistence_volume(
            &domain,
            outlet(0),
            &mixture,
            1,
            data,
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let general = field
        .advance_thermal_volume(
            &domain,
            outlet(0),
            &mixture,
            1,
            data,
            thermal_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(general.volume, original);
    assert_eq!(general.source_iterations, 1);
    assert_eq!(general.mixing_volume_source_m3_s, [0.0]);
}

#[test]
fn sensible_caloric_datum_and_relaxation_retain_the_same_physical_solution() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mut results = Vec::new();
    for (offset, relaxation) in [(0.0, 1.0), (65536.0, 0.5)] {
        let curve = sensible_chart(offset);
        let mass = [1.5];
        let states = [Some(
            curve.state_at_specific_enthalpy(offset + 2.0).unwrap(),
        )];
        let increments = [-1.5];
        let mut data = input(&mass, &states, &increments);
        data.boundary_inflow[1] = Some(curve.state_at_specific_enthalpy(offset + 3.0).unwrap());
        let mut ctl = thermal_controls();
        ctl.source_relaxation = relaxation;
        let result = field
            .advance_thermal_volume(
                &domain,
                outlet(0),
                &curve,
                0,
                data,
                ctl,
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        close(result.volume.report.requested_enthalpy_increment_j, -1.5);
        assert!(
            result.volume.report.absolute_enthalpy_defect_j
                <= ctl.volume.transport.maximum_enthalpy_defect_j
        );
        results.push(result);
    }
    let a = &results[0];
    let b = &results[1];
    assert!((a.volume.transport.mass_kg[0] - b.volume.transport.mass_kg[0]).abs() < 5e-10);
    assert!(
        (a.volume.transport.specific_enthalpy_j_kg[0]
            - (b.volume.transport.specific_enthalpy_j_kg[0] - 65536.0))
            .abs()
            < 5e-10
    );
    assert!((a.volume.volume_source_m3_s[0] - b.volume.volume_source_m3_s[0]).abs() < 5e-10);
    assert!(
        (a.volume.report.required_source_pressure_work_j
            - b.volume.report.required_source_pressure_work_j)
            .abs()
            < 5e-10
    );
    assert!(b.source_iterations > a.source_iterations);
}

#[test]
fn sensible_iteration_budget_invalid_controls_and_sealed_expansion_refuse_atomically() {
    let curve = sensible_chart(0.0);
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [1.5];
    let states = [Some(curve.state_at_specific_enthalpy(2.0).unwrap())];
    let increments = [-1.5];
    let mut data = input(&mass, &states, &increments);
    data.boundary_inflow[1] = Some(curve.state_at_specific_enthalpy(3.0).unwrap());
    let gate = CancelGate::new_clock_free();
    let run = |ctl| field.advance_thermal_volume(&domain, outlet(0), &curve, 0, data, ctl, &gate);
    for count in [0, 1] {
        let mut ctl = thermal_controls();
        ctl.max_volume_iterations = count;
        assert!(run(ctl).is_err());
    }
    for relaxation in [0.0, -1.0, 1.01, f64::NAN] {
        let mut ctl = thermal_controls();
        ctl.source_relaxation = relaxation;
        assert!(run(ctl).is_err());
    }
    assert!(
        field
            .advance_thermal_volume(
                &domain,
                [FlowFace::Wall; 6],
                &curve,
                0,
                input(&mass, &states, &[1.5]),
                thermal_controls(),
                &gate
            )
            .is_err()
    );
    assert!(
        field
            .advance_thermal_volume(
                &domain,
                outlet(0),
                &curve,
                0,
                input(&mass, &states, &[4.5]),
                thermal_controls(),
                &gate
            )
            .is_err()
    );
    let accepted = run(thermal_controls()).unwrap();
    assert_eq!(run(thermal_controls()).unwrap(), accepted);
    gate.request();
    assert!(run(thermal_controls()).is_err());
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [1.5]);
    assert_eq!(states[0].unwrap().specific_enthalpy_j_kg(), 2.0);
}

#[test]
fn sensible_zero_input_preserves_admitted_roundoff_without_a_volume_source() {
    let curve = sensible_chart(0.0);
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let mass = [f64::from_bits(1.5_f64.to_bits() + 1)];
    let states = [Some(curve.state_at_specific_enthalpy(2.0).unwrap())];
    let field = FlowField::quiescent(&domain);
    let result = field
        .advance_thermal_volume(
            &domain,
            [FlowFace::Wall; 6],
            &curve,
            0,
            input(&mass, &states, &[0.0]),
            thermal_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(result.volume.transport.mass_kg, mass);
    assert_eq!(result.volume.phase_states, states);
    assert_eq!(result.volume.transport.projection.field, field);
    assert_eq!(result.volume.volume_source_m3_s, [0.0]);
    assert_eq!(result.mixing_volume_source_m3_s, [0.0]);
    assert_eq!(result.source_iterations, 1);
}
