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
