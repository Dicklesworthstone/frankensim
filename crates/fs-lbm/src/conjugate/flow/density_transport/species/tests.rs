//! Independent material mixing, extensive conservation and atomic refusal.
use super::*;
use crate::conjugate::krylov::{bicgstab_ilu0, bicgstab_ilu0_prepared};
use crate::conjugate::{
    AdvectedMaterial, DensityProjectionControls, FlowFace, FlowField, MaterialAdvectionControls,
    MaterialAdvectionInput, MaterialVolumeSource, VoxelDomain,
};
use fs_sparse::precond::ilu0;

fn carrier_controls(dt: f64) -> MaterialAdvectionControls {
    MaterialAdvectionControls {
        projection: DensityProjectionControls {
            duration_s: dt,
            relative_tolerance: 1e-12,
            maximum_volume_defect_m3_s: 1e-10,
            maximum_energy_defect_j: 1e-8,
        },
        max_cells: 100,
        max_connections: 1000,
        max_iterations: 1000,
        relative_tolerance: 1e-12,
        maximum_mass_defect_kg: 1e-9,
        maximum_enthalpy_defect_j: 1e-7,
    }
}

fn controls<const N: usize>() -> MaterialSpeciesControls<N> {
    MaterialSpeciesControls {
        max_iterations: 1000,
        relative_tolerance: 1e-13,
        maximum_mass_defect_kg: [1e-10; N],
        maximum_partition_defect_kg: 1e-10,
        maximum_storage_bytes: 1 << 20,
    }
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(f64::MIN_POSITIVE),
        "{actual:.17e} != {expected:.17e}",
    );
}

fn carrier_input<'a>(
    mass: &'a [f64],
    h: &'a [f64],
    sources: &'a [MaterialVolumeSource],
) -> MaterialAdvectionInput<'a> {
    MaterialAdvectionInput {
        mass_kg: mass,
        specific_enthalpy_j_kg: h,
        volume_sources: sources,
        boundary_inflow: [None; 6],
    }
}

fn species_input<'a, const N: usize>(
    masses: &'a [[f64; N]],
    sources: &'a [Option<[f64; N]>],
) -> MaterialSpeciesInput<'a, N> {
    MaterialSpeciesInput {
        species_mass_kg: masses,
        source_mass_fractions: sources,
        boundary_mass_fractions: [None; 6],
    }
}

fn open_step(axis: usize, reverse: bool) -> MaterialAdvectionStep {
    let mut widths = [vec![2.0], vec![3.0], vec![2.0]];
    widths[axis] = vec![1.0, 2.0];
    widths[(axis + 1) % 3] = vec![2.0];
    widths[(axis + 2) % 3] = vec![3.0];
    let domain = VoxelDomain::graded(widths).unwrap();
    let mut faces = [FlowFace::Wall; 6];
    faces[2 * axis] = FlowFace::Fixed;
    faces[2 * axis + 1] = FlowFace::Fixed;
    let field = FlowField::from_face_velocity(&domain, faces, |a, _| {
        if a != axis {
            0.0
        } else if reverse {
            -1.0
        } else {
            1.0
        }
    })
    .unwrap();
    let mass = [12.0, 96.0];
    let h = [20.0, 80.0];
    let sources = [MaterialVolumeSource::None; 2];
    let mut input = carrier_input(&mass, &h, &sources);
    input.boundary_inflow[2 * axis + usize::from(reverse)] = Some(AdvectedMaterial {
        density_kg_m3: 10.0,
        specific_enthalpy_j_kg: 100.0,
    });
    field
        .project_and_advect_material(
            &domain,
            faces,
            input,
            carrier_controls(0.5),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
}

#[test]
fn graded_density_mixing_matches_partial_mass_balances_on_each_axis_and_backflow() {
    // V=[6,12], dt*Q=3. Solve independent partial-density storage equations:
    // (V+dtQ)*rho_k = initial species mass + incoming species mass.
    for axis in 0..3 {
        for reverse in [false, true] {
            let carrier = open_step(axis, reverse);
            let masses = [[9.0, 3.0, 0.0], [24.0, 72.0, 0.0]];
            let inlet = [0.2, 0.8, 0.0];
            let mut input = species_input(&masses, &[None; 2]);
            input.boundary_mass_fractions[2 * axis + usize::from(reverse)] = Some(inlet);
            let next = carrier
                .advect_species(input, controls(), &CancelGate::new_clock_free())
                .unwrap();
            for k in 0..2 {
                let expected = if reverse {
                    let right = (masses[1][k] + 30.0 * inlet[k]) / 15.0;
                    [(masses[0][k] + 3.0 * right) * (6.0 / 9.0), right * 12.0]
                } else {
                    let left = (masses[0][k] + 30.0 * inlet[k]) / 9.0;
                    [left * 6.0, (masses[1][k] + 3.0 * left) * (12.0 / 15.0)]
                };
                for (c, expected) in expected.into_iter().enumerate() {
                    close(next.species_mass_kg[c][k], expected, 3e-13);
                }
                close(
                    next.report.boundary_incoming_mass_kg[k],
                    30.0 * inlet[k],
                    3e-13,
                );
                assert_eq!(next.report.source_incoming_mass_kg[k], 0.0);
            }
            assert_eq!(next.species_mass_kg[0][2], 0.0);
            assert_eq!(next.species_mass_kg[1][2], 0.0);
            assert_eq!(next.report.iterations[2], 0);
            for (c, values) in next.species_mass_kg.iter().enumerate() {
                close(values.iter().sum(), carrier.mass_kg[c], 3e-13);
            }
        }
    }
}

#[test]
fn coincident_boundary_and_local_inlets_keep_their_distinct_compositions() {
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
        FlowField::from_face_velocity(&domain, faces, |a, _| if a == 0 { 1.0 } else { 0.0 })
            .unwrap();
    let sources = [MaterialVolumeSource::Inflow {
        volume_rate_m3_s: 0.5,
        material: AdvectedMaterial {
            density_kg_m3: 8.0,
            specific_enthalpy_j_kg: 20.0,
        },
    }];
    let mut input = carrier_input(&[2.0], &[20.0], &sources);
    input.boundary_inflow[0] = Some(AdvectedMaterial {
        density_kg_m3: 4.0,
        specific_enthalpy_j_kg: 20.0,
    });
    let carrier = field
        .project_and_advect_material(
            &domain,
            faces,
            input,
            carrier_controls(0.5),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // V=1, dt*Q_out=3/4. Each distinct stream supplies 2 kg, so every
    // endpoint inventory is 2/(1+3/4)=8/7 kg; no 1/3 composition is authored.
    let masses = [[0.0, 0.0, 2.0]];
    let origins = [Some([1.0, 0.0, 0.0])];
    let mut input = species_input(&masses, &origins);
    input.boundary_mass_fractions[0] = Some([0.0, 1.0, 0.0]);
    let next = carrier
        .advect_species(input, controls(), &CancelGate::new_clock_free())
        .unwrap();
    close(carrier.mass_kg[0], 24.0 / 7.0, 3e-13);
    for k in 0..3 {
        close(next.species_mass_kg[0][k], 8.0 / 7.0, 3e-13);
        close(next.report.outgoing_mass_kg[k], 6.0 / 7.0, 3e-13);
    }
    assert_eq!(next.report.source_incoming_mass_kg, [2.0, 0.0, 0.0]);
    assert_eq!(next.report.boundary_incoming_mass_kg, [0.0, 2.0, 0.0]);
    let without_source = species_input(&masses, &[None]);
    assert!(
        carrier
            .advect_species(
                MaterialSpeciesInput {
                    boundary_mass_fractions: input.boundary_mass_fractions,
                    ..without_source
                },
                controls(),
                &CancelGate::new_clock_free(),
            )
            .is_err()
    );
    let mut without_boundary = input;
    without_boundary.boundary_mass_fractions = [None; 6];
    assert!(
        carrier
            .advect_species(without_boundary, controls(), &CancelGate::new_clock_free())
            .is_err()
    );
}

#[test]
fn trace_mass_survives_native_solve_without_seeding_absent_species() {
    let carrier = open_step(0, false);
    let trace = 1e-200;
    let masses = [[12.0, 12.0 * trace, 0.0], [96.0, 192.0 * trace, 0.0]];
    let mut input = species_input(&masses, &[None; 2]);
    input.boundary_mass_fractions[0] = Some([1.0, 3.0 * trace, 0.0]);
    let mut ctl = controls();
    ctl.maximum_mass_defect_kg[1] = 1e-210;
    let next = carrier
        .advect_species(input, ctl, &CancelGate::new_clock_free())
        .unwrap();
    let left = (12.0 * trace + 90.0 * trace) / 9.0;
    let right = (192.0 * trace + 3.0 * left) / 15.0;
    close(next.species_mass_kg[0][1], 6.0 * left, 5e-13);
    close(next.species_mass_kg[1][1], 12.0 * right, 5e-13);
    assert_eq!(next.species_mass_kg[0][2], 0.0);
    assert_eq!(next.species_mass_kg[1][2], 0.0);
    assert_eq!(next.report.absolute_mass_defect_kg[2], 0.0);
    assert!(next.report.absolute_mass_defect_kg[1] < ctl.maximum_mass_defect_kg[1]);
}

#[test]
fn pure_component_on_nonbinary_volumes_retains_original_extensive_inventory() {
    let domain = VoxelDomain::graded([vec![0.1, 0.2, 0.3], vec![0.7], vec![0.9]]).unwrap();
    let field = FlowField::quiescent(&domain);
    let masses = [0.13, 1.7, 0.29];
    let carrier = field
        .project_and_advect_material(
            &domain,
            [FlowFace::Wall; 6],
            carrier_input(&masses, &[20.0; 3], &[MaterialVolumeSource::None; 3]),
            carrier_controls(0.4),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(carrier.mass_kg, masses);
    let species = masses.map(|mass| [0.0, mass, 0.0]);
    let next = carrier
        .advect_species(
            species_input(&species, &[None; 3]),
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(next.species_mass_kg, species);
    assert_eq!(next.report.endpoint_partition_defect_kg, 0.0);
    assert_eq!(next.report.absolute_mass_defect_kg, [0.0; 3]);
    // The represented non-binary volume does not authorize replacing m by
    // an independently rounded rho*V or repairing a pure fraction afterward.
    assert_eq!(next.report.initial_partition_defect_kg, 0.0);
}

#[test]
fn flowing_pure_species_admits_carrier_roundoff_without_normalizing_fractions() {
    let domain = VoxelDomain::graded([vec![0.1, 0.2], vec![0.7], vec![0.9]]).unwrap();
    let faces = [
        FlowFace::Fixed,
        FlowFace::Fixed,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let field =
        FlowField::from_face_velocity(&domain, faces, |a, _| if a == 0 { 0.4 } else { 0.0 })
            .unwrap();
    let masses = [0.13, 1.7];
    let mut input = carrier_input(&masses, &[20.0; 2], &[MaterialVolumeSource::None; 2]);
    input.boundary_inflow[0] = Some(AdvectedMaterial {
        density_kg_m3: 1.3,
        specific_enthalpy_j_kg: 20.0,
    });
    let carrier = field
        .project_and_advect_material(
            &domain,
            faces,
            input,
            carrier_controls(0.17),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let pure = masses.map(|mass| [0.0, mass, 0.0]);
    let mut input = species_input(&pure, &[None; 2]);
    input.boundary_mass_fractions[0] = Some([0.0, 1.0, 0.0]);
    let next = carrier
        .advect_species(input, controls(), &CancelGate::new_clock_free())
        .unwrap();
    let mut actual_partition = 0.0;
    for c in 0..2 {
        let returned = next.species_mass_kg[c];
        assert_eq!(returned[0], 0.0);
        assert_eq!(returned[2], 0.0);
        close(returned[1], carrier.mass_kg[c], 5e-13);
        actual_partition += (returned[1] - carrier.mass_kg[c]).abs();
        // A derived pure fraction may leave exact [0,1] by an ulp. The
        // declared extensive allowance, not a clamp, decides acceptance.
        close(returned[1] / carrier.mass_kg[c], 1.0, 5e-13);
    }
    assert_eq!(next.report.endpoint_partition_defect_kg, actual_partition);
    assert_eq!(next.report.absolute_mass_defect_kg[0], 0.0);
    assert_eq!(next.report.absolute_mass_defect_kg[2], 0.0);
}

fn circulation(domain: &VoxelDomain) -> FlowField {
    FlowField::from_face_velocity(domain, [FlowFace::Wall; 6], |axis, p| match axis {
        0 => {
            if p[1] < 1.0 {
                1.0
            } else {
                -1.0
            }
        }
        1 => {
            if p[0] < 1.0 {
                -1.0
            } else {
                1.0
            }
        }
        _ => 0.0,
    })
    .unwrap()
}

#[test]
fn repeated_closed_circulation_conserves_each_species_and_permuting_species_commutes() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mut field = circulation(&domain);
    let mut mass = vec![1.0, 2.0, 3.0, 4.0];
    let mut enthalpy = vec![20.0; 4];
    let mut species = vec![
        [1.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.0, 0.0, 3.0],
        [1.0, 1.0, 2.0],
    ];
    for _ in 0..16 {
        let carrier = field
            .project_and_advect_material(
                &domain,
                [FlowFace::Wall; 6],
                carrier_input(&mass, &enthalpy, &[MaterialVolumeSource::None; 4]),
                carrier_controls(0.25),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        let next = carrier
            .advect_species(
                species_input(&species, &[None; 4]),
                controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        let reordered: Vec<_> = species.iter().map(|v| [v[2], v[0], v[1]]).collect();
        let permuted = carrier
            .advect_species(
                species_input(&reordered, &[None; 4]),
                controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        // Independent cyclic partial-density equations: x_i = 0.8*m_i
        // + 0.2*x_previous. Unroll the four-cell ring analytically.
        let ring = [0, 1, 3, 2];
        for (position, &c) in ring.iter().enumerate() {
            let values = &next.species_mass_kg[c];
            assert_eq!(
                permuted.species_mass_kg[c],
                [values[2], values[0], values[1]]
            );
            close(values.iter().sum(), carrier.mass_kg[c], 5e-13);
            for k in 0..3 {
                let mut sum = 0.0;
                let mut weight = 1.0;
                for back in 0..4 {
                    sum += weight * species[ring[(position + 4 - back) % 4]][k];
                    weight *= 0.2;
                }
                close(values[k], 0.8 * sum / (1.0 - 0.0016), 5e-13);
            }
        }
        for (k, total) in [2.0, 3.0, 5.0].into_iter().enumerate() {
            close(next.report.mass_after_kg[k], total, 5e-13);
            assert_eq!(next.report.boundary_incoming_mass_kg[k], 0.0);
            assert_eq!(next.report.source_incoming_mass_kg[k], 0.0);
            assert_eq!(next.report.outgoing_mass_kg[k], 0.0);
        }
        field = carrier.projection.field;
        mass = carrier.mass_kg;
        enthalpy = carrier.specific_enthalpy_j_kg;
        species = next.species_mass_kg;
    }
    assert!(species.iter().all(|v| v.iter().all(|&m| m > 0.0)));
}

#[test]
fn shared_factor_keeps_original_native_solver_arithmetic_and_zero_rhs_path() {
    let carrier = open_step(1, false);
    let a = &carrier.mass_transport.operator;
    let factor = ilu0(a).unwrap();
    for b in [[0.2, 0.8], [0.7, 0.1], [0.0, 0.0]] {
        let mut native = [0.31, 0.17];
        let mut shared = native;
        let gate = CancelGate::new_clock_free();
        let expected = bicgstab_ilu0("parity", a, &b, &mut native, 1e-13, 100, &gate).unwrap();
        let mut constructed = false;
        let actual =
            bicgstab_ilu0_prepared("parity", a, &b, &mut shared, 1e-13, 100, &gate, || {
                constructed = true;
                Ok(&factor)
            })
            .unwrap();
        assert_eq!(native, shared);
        assert_eq!(expected.iterations, actual.iterations);
        assert_eq!(expected.relative_residual, actual.relative_residual);
        assert_eq!(constructed, b != [0.0, 0.0]);
    }
}

#[test]
fn late_refusal_cancellation_and_retry_never_publish_a_partial_constituent_state() {
    let carrier = open_step(0, false);
    let before = carrier.clone();
    let masses = [[9.0, 3.0], [24.0, 72.0]];
    let mut input = species_input(&masses, &[None; 2]);
    input.boundary_mass_fractions[0] = Some([0.2, 0.8]);
    let original_input = masses;
    let gate = CancelGate::new_clock_free();
    let mut completed = 0;
    let cancelled = carrier.advect_species_checked(input, controls(), &gate, |count| {
        completed = count;
        if count == 1 {
            gate.request();
        }
        poll(&gate)
    });
    assert_eq!(completed, 1);
    assert_eq!(cancelled, Err(ChtError::Cancelled));
    assert_eq!(carrier, before);
    assert_eq!(masses, original_input);

    // Refuse a real unrepresentable positive species only AFTER one complete
    // species has run; no seeded epsilon or normalized bulk replacement.
    let tiny = 12.0 * f64::from_bits(1);
    let late_masses = [[12.0, tiny], [96.0, 0.0]];
    let mut late = species_input(&late_masses, &[None; 2]);
    late.boundary_mass_fractions[0] = Some([1.0, 0.0]);
    let gate = CancelGate::new_clock_free();
    let mut completed = 0;
    let refused = carrier.advect_species_checked(late, controls(), &gate, |count| {
        completed = count;
        poll(&gate)
    });
    assert_eq!(completed, 1);
    assert!(matches!(
        refused,
        Err(ChtError::InvalidInput {
            field: "flow.material_species",
            ..
        })
    ));
    assert_eq!(carrier, before);
    assert_eq!(late_masses[0][1].to_bits(), tiny.to_bits());

    let mut too_small = controls();
    too_small.maximum_storage_bytes = 0;
    assert!(
        carrier
            .advect_species(input, too_small, &CancelGate::new_clock_free())
            .is_err()
    );
    let bad_partition = [[8.0, 3.0], [24.0, 72.0]];
    assert!(
        carrier
            .advect_species(
                MaterialSpeciesInput {
                    species_mass_kg: &bad_partition,
                    ..input
                },
                controls(),
                &CancelGate::new_clock_free(),
            )
            .is_err()
    );
    let mut changed_carrier = carrier.clone();
    changed_carrier.mass_kg[0] += 1.0;
    assert!(
        changed_carrier
            .advect_species(input, controls(), &CancelGate::new_clock_free())
            .is_err()
    );

    let retried = carrier
        .advect_species(input, controls(), &CancelGate::new_clock_free())
        .unwrap();
    let fresh = open_step(0, false)
        .advect_species(input, controls(), &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(retried, fresh);
}
