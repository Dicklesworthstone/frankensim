//! Independent source/energy checks through the original material and species owners.
use super::*;
use fs_lbm::conjugate::{
    GasSpeciesAdvectionInput, GasSpeciesVolumeSource, MaterialSpeciesControls,
    MaterialSpeciesInput, ThermochemicalSpeciesAdvectionControls,
};

fn advection_controls(dt: f64) -> ThermochemicalSpeciesAdvectionControls<2> {
    ThermochemicalSpeciesAdvectionControls {
        material: MaterialAdvectionControls {
            projection: DensityProjectionControls {
                duration_s: dt,
                relative_tolerance: 1e-12,
                maximum_volume_defect_m3_s: 1e-10,
                maximum_energy_defect_j: 1e-8,
            },
            max_cells: 16,
            max_connections: 128,
            max_iterations: 256,
            relative_tolerance: 1e-12,
            maximum_mass_defect_kg: 1e-10,
            maximum_enthalpy_defect_j: 1e-6,
        },
        species: MaterialSpeciesControls {
            max_iterations: 256,
            relative_tolerance: 1e-13,
            maximum_mass_defect_kg: [1e-10; 2],
            maximum_partition_defect_kg: 1e-10,
            maximum_storage_bytes: 1 << 20,
        },
        maximum_source_evaluations: 32,
        maximum_additional_storage_bytes: 1 << 20,
        maximum_enthalpy_defect_j: ENERGY_BUDGET,
    }
}

#[test]
fn g1_source_bound_boundary_and_local_inflow_reconstruct_analytic_mixture_temperature() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let faces = [
        FlowFace::Fixed,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let flow =
        FlowField::from_face_velocity(&domain, faces, |axis, _| if axis == 0 { 0.4 } else { 0.0 })
            .unwrap();
    let original_flow = flow.clone();
    let names = ids();
    let context = source::point(400.0, 100_000.0);
    let gate = CancelGate::new_clock_free();
    let diffusivity = [0.0];
    let mass = [2.0];
    let constituents = [[1.5, 0.5]];
    let mut previous_temperature: Option<f64> = None;
    for a7 in [-4000.0, 0.0, 10000.0] {
        let bank = source::moist_model(0.5, a7);
        let gas =
            GasSpeciesEnthalpySource::new(&bank, &names, &context, &diffusivity, energy(), &gate)
                .unwrap();
        let incoming = gas.advection_inlet([1.0, 0.0], 600.0, &gate).unwrap();
        let injection = gas.advection_inlet([0.0, 1.0], 450.0, &gate).unwrap();
        assert_eq!(incoming.mass_fractions(), &[1.0, 0.0]);
        assert_eq!(injection.mass_fractions(), &[0.0, 1.0]);
        assert_eq!(incoming.temperature_k(), 600.0);
        let source_h = partial_h(320.0, a7);
        let h = [(1.5 * source_h[0] + 0.5 * source_h[1]) / 2.0];
        let sources = [GasSpeciesVolumeSource::Inflow {
            volume_rate_m3_s: 0.2,
            inlet: injection,
        }];
        let input = GasSpeciesAdvectionInput {
            mass_kg: &mass,
            species_mass_kg: &constituents,
            specific_enthalpy_j_kg: &h,
            volume_sources: &sources,
            boundary_inflow: [Some(incoming), None, None, None, None, None],
        };
        let controls = advection_controls(0.25);
        let result = flow
            .project_and_advect_thermochemical_species(&domain, faces, input, &gas, controls, &gate)
            .unwrap();

        // The original volume projection gives Qout=.4+.2=.6. Its backward-
        // Euler donor equation is (1+dt*Qout/V) m_k,new = m_k,old+dt*Qin*rho*Y.
        let incoming_mass = 0.25 * 0.4 * 100_000.0 * MOLAR_MASS[0] / (R * 600.0);
        let injected_mass = 0.25 * 0.2 * 100_000.0 * MOLAR_MASS[1] / (R * 450.0);
        let denominator = 1.0 + 0.25 * 0.6;
        let expected_species = [
            (1.5 + incoming_mass) / denominator,
            (0.5 + injected_mass) / denominator,
        ];
        let expected_mass = expected_species.iter().sum::<f64>();
        let total_h = (2.0 * h[0]
            + incoming_mass * partial_h(600.0, a7)[0]
            + injected_mass * partial_h(450.0, a7)[1])
            / denominator;
        let expected_t = temperature(total_h, expected_species, a7);
        for k in 0..2 {
            near(
                result.species.species_mass_kg[0][k],
                expected_species[k],
                1e-12,
            );
        }
        near(result.carrier.mass_kg[0], expected_mass, 1e-12);
        near(
            result.carrier.specific_enthalpy_j_kg[0],
            total_h / expected_mass,
            1e-6,
        );
        near(result.state.temperatures_k()[0], expected_t, 1e-8);
        near(
            result.species.report.boundary_incoming_mass_kg[0],
            incoming_mass,
            1e-12,
        );
        near(
            result.species.report.source_incoming_mass_kg[1],
            injected_mass,
            1e-12,
        );
        assert_eq!(result.species.report.boundary_incoming_mass_kg[1], 0.0);
        assert_eq!(result.species.report.source_incoming_mass_kg[0], 0.0);
        assert_eq!(result.report.source_evaluations, 2);
        assert!(result.report.absolute_enthalpy_defect_j <= ENERGY_BUDGET);
        assert_eq!(result.state.species(), &names);
        assert_eq!(result.state.source().source_packs(), bank.source_packs());
        if let Some(previous) = previous_temperature {
            near(result.state.temperatures_k()[0], previous, 1e-8);
        }
        previous_temperature = Some(result.state.temperatures_k()[0]);

        // Exact native parity: the adapter must not reconstruct different mass
        // routes, replace donor fractions, or apply another thermal deposit.
        let native_sources = [MaterialVolumeSource::Inflow {
            volume_rate_m3_s: 0.2,
            material: injection.advected_material(),
        }];
        let native = flow
            .project_and_advect_material(
                &domain,
                faces,
                MaterialAdvectionInput {
                    mass_kg: &mass,
                    specific_enthalpy_j_kg: &h,
                    volume_sources: &native_sources,
                    boundary_inflow: [
                        Some(incoming.advected_material()),
                        None,
                        None,
                        None,
                        None,
                        None,
                    ],
                },
                controls.material,
                &gate,
            )
            .unwrap();
        let native_species = native
            .advect_species(
                MaterialSpeciesInput {
                    species_mass_kg: &constituents,
                    source_mass_fractions: &[Some([0.0, 1.0])],
                    boundary_mass_fractions: [Some([1.0, 0.0]), None, None, None, None, None],
                },
                controls.species,
                &gate,
            )
            .unwrap();
        assert_eq!(result.carrier, native);
        assert_eq!(result.species, native_species);
    }
    assert_eq!(flow, original_flow);
    assert_eq!(mass, [2.0]);
    assert_eq!(constituents, [[1.5, 0.5]]);
}

#[test]
fn g1_reverse_flow_uses_its_declared_positive_side_inlet_and_supports_continuation() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let faces = [
        FlowFace::Fixed,
        FlowFace::Fixed,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let flow = FlowField::from_face_velocity(
        &domain,
        faces,
        |axis, _| {
            if axis == 0 { -0.25 } else { 0.0 }
        },
    )
    .unwrap();
    let bank = source::moist_model(0.5, -4000.0);
    let names = ids();
    let context = source::point(400.0, 100_000.0);
    let diffusivity = [0.0];
    let gate = CancelGate::new_clock_free();
    let gas = GasSpeciesEnthalpySource::new(&bank, &names, &context, &diffusivity, energy(), &gate)
        .unwrap();
    let inlet = gas.advection_inlet([0.0, 1.0], 500.0, &gate).unwrap();
    let mass = [1.0];
    let constituents = [[1.0, 0.0]];
    let h = [partial_h(350.0, -4000.0)[0]];
    let sources = [GasSpeciesVolumeSource::None];
    let input = GasSpeciesAdvectionInput {
        mass_kg: &mass,
        species_mass_kg: &constituents,
        specific_enthalpy_j_kg: &h,
        volume_sources: &sources,
        boundary_inflow: [None, Some(inlet), None, None, None, None],
    };
    let result = flow
        .project_and_advect_thermochemical_species(
            &domain,
            faces,
            input,
            &gas,
            advection_controls(0.2),
            &gate,
        )
        .unwrap();
    let dm = 0.05 * 100_000.0 * MOLAR_MASS[1] / (R * 500.0);
    let expected = [1.0 / 1.05, dm / 1.05];
    let total_h = (h[0] + dm * partial_h(500.0, -4000.0)[1]) / 1.05;
    for k in 0..2 {
        near(result.species.species_mass_kg[0][k], expected[k], 1e-12);
    }
    near(
        result.state.temperatures_k()[0],
        temperature(total_h, expected, -4000.0),
        1e-8,
    );
    assert_eq!(result.species.report.boundary_incoming_mass_kg[0], 0.0);
    near(
        result.species.report.boundary_incoming_mass_kg[1],
        dm,
        1e-12,
    );
    let next = result
        .carrier
        .projection
        .field
        .project_and_advect_thermochemical_species(
            &domain,
            faces,
            GasSpeciesAdvectionInput {
                mass_kg: &result.carrier.mass_kg,
                species_mass_kg: &result.species.species_mass_kg,
                specific_enthalpy_j_kg: &result.carrier.specific_enthalpy_j_kg,
                volume_sources: &sources,
                boundary_inflow: input.boundary_inflow,
            },
            &gas,
            advection_controls(0.2),
            &gate,
        )
        .unwrap();
    let expected_next = [expected[0] / 1.05, (expected[1] + dm) / 1.05];
    for k in 0..2 {
        near(next.species.species_mass_kg[0][k], expected_next[k], 2e-12);
    }
    assert!(next.state.temperatures_k()[0] > result.state.temperatures_k()[0]);
    assert!(
        flow.project_and_advect_thermochemical_species(
            &domain,
            faces,
            GasSpeciesAdvectionInput {
                boundary_inflow: [Some(inlet), None, None, None, None, None],
                ..input
            },
            &gas,
            advection_controls(0.2),
            &gate,
        )
        .is_err()
    );
}

#[test]
fn g0_source_binding_partition_resource_and_cancellation_refusals_are_atomic() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let faces = [FlowFace::Wall; 6];
    let flow = FlowField::from_face_velocity(&domain, faces, |_, _| 0.0).unwrap();
    let before = flow.clone();
    let bank = source::moist_model(0.5, 0.0);
    let other_bank = source::moist_model(0.5, 0.0);
    let names = ids();
    let context = source::point(400.0, 100_000.0);
    let other_context = source::point(400.0, 100_000.0);
    let diffusivity = [0.0];
    let gate = CancelGate::new_clock_free();
    let gas = GasSpeciesEnthalpySource::new(&bank, &names, &context, &diffusivity, energy(), &gate)
        .unwrap();
    let foreign =
        GasSpeciesEnthalpySource::new(&other_bank, &names, &context, &diffusivity, energy(), &gate)
            .unwrap();
    let foreign_context =
        GasSpeciesEnthalpySource::new(&bank, &names, &other_context, &diffusivity, energy(), &gate)
            .unwrap();
    assert!(gas.advection_inlet([0.6, 0.6], 400.0, &gate).is_err());
    assert!(gas.advection_inlet([1.0, 0.0], 199.0, &gate).is_err());
    let mass = [1.0];
    let constituents = [[1.0, 0.0]];
    let h = [partial_h(350.0, 0.0)[0]];
    let sources = [GasSpeciesVolumeSource::None];
    let input = GasSpeciesAdvectionInput {
        mass_kg: &mass,
        species_mass_kg: &constituents,
        specific_enthalpy_j_kg: &h,
        volume_sources: &sources,
        boundary_inflow: [None; 6],
    };
    for inlet in [
        foreign.advection_inlet([1.0, 0.0], 400.0, &gate).unwrap(),
        foreign_context
            .advection_inlet([1.0, 0.0], 400.0, &gate)
            .unwrap(),
    ] {
        assert!(
            flow.project_and_advect_thermochemical_species(
                &domain,
                faces,
                GasSpeciesAdvectionInput {
                    boundary_inflow: [Some(inlet), None, None, None, None, None],
                    ..input
                },
                &gas,
                advection_controls(0.2),
                &gate,
            )
            .is_err()
        );
    }
    for choice in 0..3 {
        let mut controls = advection_controls(0.2);
        match choice {
            0 => controls.maximum_source_evaluations = 1,
            1 => controls.maximum_additional_storage_bytes = 0,
            _ => controls.material.maximum_enthalpy_defect_j = f64::NAN,
        }
        assert!(
            flow.project_and_advect_thermochemical_species(
                &domain, faces, input, &gas, controls, &gate,
            )
            .is_err()
        );
    }
    assert!(
        flow.project_and_advect_thermochemical_species(
            &domain,
            faces,
            GasSpeciesAdvectionInput {
                species_mass_kg: &[[1.0001, 0.0]],
                ..input
            },
            &gas,
            advection_controls(0.2),
            &gate,
        )
        .is_err()
    );
    assert!(
        flow.project_and_advect_thermochemical_species(
            &domain,
            faces,
            GasSpeciesAdvectionInput {
                mass_kg: &[f64::MAX],
                species_mass_kg: &[[f64::MIN_POSITIVE, f64::MAX]],
                ..input
            },
            &gas,
            advection_controls(0.2),
            &gate,
        )
        .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        flow.project_and_advect_thermochemical_species(
            &domain,
            faces,
            input,
            &gas,
            advection_controls(0.2),
            &cancelled,
        ),
        Err(ChtError::Cancelled),
    ));
    assert_eq!(flow, before);
    assert_eq!(mass, [1.0]);
    assert_eq!(constituents, [[1.0, 0.0]]);
    assert_eq!(input.specific_enthalpy_j_kg, h);
}

#[test]
fn g3_closed_quiescent_source_keeps_absent_species_and_solid_entries() {
    let mut domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    domain.set(2, 0, 0, Voxel::Solid(0));
    let faces = [FlowFace::Wall; 6];
    let flow = FlowField::from_face_velocity(&domain, faces, |_, _| 0.0).unwrap();
    let bank = source::moist_model(0.5, -4000.0);
    let names = ids();
    let context = source::point(400.0, 100_000.0);
    let diffusivity = [0.0; 3];
    let gate = CancelGate::new_clock_free();
    let gas = GasSpeciesEnthalpySource::new(&bank, &names, &context, &diffusivity, energy(), &gate)
        .unwrap();
    let masses = [1.0, 2.0, 0.0];
    let species = [[1.0, 0.0], [0.0, 2.0], [0.0; 2]];
    let h = [
        partial_h(300.0, -4000.0)[0],
        partial_h(450.0, -4000.0)[1],
        0.0,
    ];
    let sources = [GasSpeciesVolumeSource::None; 3];
    let result = flow
        .project_and_advect_thermochemical_species(
            &domain,
            faces,
            GasSpeciesAdvectionInput {
                mass_kg: &masses,
                species_mass_kg: &species,
                specific_enthalpy_j_kg: &h,
                volume_sources: &sources,
                boundary_inflow: [None; 6],
            },
            &gas,
            advection_controls(0.2),
            &gate,
        )
        .unwrap();
    assert_eq!(result.carrier.mass_kg, masses);
    assert_eq!(result.species.species_mass_kg, species);
    assert_eq!(result.carrier.specific_enthalpy_j_kg, h);
    near(result.state.temperatures_k()[0], 300.0, 1e-8);
    near(result.state.temperatures_k()[1], 450.0, 1e-8);
    assert_eq!(result.state.temperatures_k()[2], 0.0);
    assert!(result.state.transport(2).is_none());
    assert_eq!(result.report.source_evaluations, 4);
}

#[test]
fn g1_signed_dilatation_redistributes_original_gas_without_local_material_sources() {
    // Unequal volumes distinguish integrated m3/s from a per-volume strain rate.
    let domain = VoxelDomain::graded([vec![1.0, 2.0], vec![1.0], vec![1.0]]).unwrap();
    let faces = [FlowFace::Wall; 6];
    let flow = FlowField::quiescent(&domain);
    let before = flow.clone();
    let mass = [2.0, 8.0];
    let constituents = [[1.5, 0.5], [2.0, 6.0]];
    let bank = source::moist_model(0.5, -4000.0);
    let names = ids();
    let context = source::point(400.0, 100_000.0);
    let diffusivity = [0.0; 2];
    let gate = CancelGate::new_clock_free();
    let gas = GasSpeciesEnthalpySource::new(&bank, &names, &context, &diffusivity, energy(), &gate)
        .unwrap();
    let h = [
        enthalpy(&bank, constituents[0], 320.0),
        enthalpy(&bank, constituents[1], 480.0),
    ];
    let sources = [
        GasSpeciesVolumeSource::Dilatation {
            volume_rate_m3_s: 1.0,
        },
        GasSpeciesVolumeSource::Dilatation {
            volume_rate_m3_s: -1.0,
        },
    ];
    let input = GasSpeciesAdvectionInput {
        mass_kg: &mass,
        species_mass_kg: &constituents,
        specific_enthalpy_j_kg: &h,
        volume_sources: &sources,
        boundary_inflow: [None; 6],
    };
    let mut controls = advection_controls(0.5);
    controls.material.max_connections = 1;
    controls.maximum_source_evaluations = 4;
    let result = flow
        .project_and_advect_thermochemical_species(&domain, faces, input, &gas, controls, &gate)
        .unwrap();

    // q=dt*Q=.5 m3 moves donor material. No local rho*Q source or sink exists:
    // M0_new=M0/(1+q/V0), transfer=q*M0_new/V0, M1_new=M1+transfer.
    let q = 0.5;
    let donor_mass = mass[0] / (1.0 + q / domain.volume(0));
    let transferred = q * donor_mass / domain.volume(0);
    let expected_mass = [donor_mass, mass[1] + transferred];
    let expected_species = [
        constituents[0].map(|m| m / mass[0] * donor_mass),
        core::array::from_fn(|k| constituents[1][k] + transferred * constituents[0][k] / mass[0]),
    ];
    let expected_h = [
        h[0],
        (mass[1] * h[1] + transferred * h[0]) / expected_mass[1],
    ];
    near(result.carrier.projection.field.flux_x(1, 0, 0), 1.0, 1e-12);
    assert_eq!(result.carrier.report.mass_routes, 1);
    for c in 0..2 {
        near(result.carrier.mass_kg[c], expected_mass[c], 1e-11);
        near(
            result.carrier.specific_enthalpy_j_kg[c],
            expected_h[c],
            1e-6,
        );
        for k in 0..2 {
            near(
                result.species.species_mass_kg[c][k],
                expected_species[c][k],
                1e-11,
            );
        }
        near(
            result.state.temperatures_k()[c],
            temperature(
                expected_mass[c] * expected_h[c],
                expected_species[c],
                -4000.0,
            ),
            1e-8,
        );
    }
    near(
        result.carrier.mass_kg.iter().sum(),
        mass.iter().sum(),
        1e-11,
    );
    near(
        result
            .carrier
            .mass_kg
            .iter()
            .zip(&result.carrier.specific_enthalpy_j_kg)
            .map(|(m, h)| m * h)
            .sum(),
        mass.iter().zip(h).map(|(m, h)| m * h).sum(),
        ENERGY_BUDGET,
    );
    assert_eq!(result.carrier.report.incoming_mass_kg, 0.0);
    assert_eq!(result.carrier.report.outgoing_mass_kg, 0.0);
    near(result.carrier.report.enthalpy_change_j, 0.0, 1e-6);
    assert_eq!(result.carrier.report.carried_enthalpy_j, 0.0);
    assert_eq!(result.species.report.source_incoming_mass_kg, [0.0; 2]);
    assert_eq!(result.species.report.boundary_incoming_mass_kg, [0.0; 2]);
    assert_eq!(result.species.report.outgoing_mass_kg, [0.0; 2]);
    for k in 0..2 {
        near(
            result.species.report.mass_after_kg[k],
            result.species.report.mass_before_kg[k],
            1e-11,
        );
    }
    assert_eq!(result.report.source_evaluations, 4);
    assert!(result.report.absolute_enthalpy_defect_j <= ENERGY_BUDGET);

    // A sealed component cannot silently remove the incompatible mean volume.
    for rate in [1.0, -1.0] {
        let incompatible = [
            GasSpeciesVolumeSource::Dilatation {
                volume_rate_m3_s: rate,
            },
            GasSpeciesVolumeSource::None,
        ];
        assert!(
            flow.project_and_advect_thermochemical_species(
                &domain,
                faces,
                GasSpeciesAdvectionInput {
                    volume_sources: &incompatible,
                    ..input
                },
                &gas,
                controls,
                &gate,
            )
            .is_err()
        );
    }
    assert_eq!(flow, before);
    assert_eq!(mass, [2.0, 8.0]);
    assert_eq!(constituents, [[1.5, 0.5], [2.0, 6.0]]);
    assert_eq!(input.specific_enthalpy_j_kg, h);
}

#[test]
fn g1_open_dilatation_moves_only_real_boundary_material_in_both_directions() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let mut faces = [FlowFace::Wall; 6];
    faces[1] = FlowFace::Free;
    let flow = FlowField::quiescent(&domain);
    let before = flow.clone();
    let mass = [2.0];
    let constituents = [[1.5, 0.5]];
    let bank = source::moist_model(0.5, -4000.0);
    let names = ids();
    let context = source::point(400.0, 100_000.0);
    let diffusivity = [0.0];
    let gate = CancelGate::new_clock_free();
    let gas = GasSpeciesEnthalpySource::new(&bank, &names, &context, &diffusivity, energy(), &gate)
        .unwrap();
    let h = [enthalpy(&bank, constituents[0], 350.0)];
    for rate in [0.4_f64, -0.4] {
        let sources = [GasSpeciesVolumeSource::Dilatation {
            volume_rate_m3_s: rate,
        }];
        let mut input = GasSpeciesAdvectionInput {
            mass_kg: &mass,
            species_mass_kg: &constituents,
            specific_enthalpy_j_kg: &h,
            volume_sources: &sources,
            boundary_inflow: [None; 6],
        };
        let mut controls = advection_controls(0.25);
        controls.material.max_connections = 1;
        controls.maximum_source_evaluations = 2;
        if rate < 0.0 {
            assert!(matches!(
                flow.project_and_advect_thermochemical_species(
                    &domain, faces, input, &gas, controls, &gate,
                ),
                Err(ChtError::InvalidInput {
                    field: "flow.material_advection",
                    ..
                })
            ));
            input.boundary_inflow[1] = Some(gas.advection_inlet([0.0, 1.0], 500.0, &gate).unwrap());
        }
        let result = flow
            .project_and_advect_thermochemical_species(&domain, faces, input, &gas, controls, &gate)
            .unwrap();
        let q = 0.25 * rate.abs();
        let (expected_species, total_h, incoming_mass, outgoing_mass) = if rate > 0.0 {
            // Expansion exports endpoint donor density and h, with no injection.
            let denominator = 1.0 + q / domain.volume(0);
            let species = constituents[0].map(|m| m / denominator);
            let remaining = mass[0] / denominator;
            (species, remaining * h[0], 0.0, mass[0] - remaining)
        } else {
            // Contraction draws the explicitly authored source-bound water inlet.
            let incoming = q * 100_000.0 * MOLAR_MASS[1] / (R * 500.0);
            (
                [constituents[0][0], constituents[0][1] + incoming],
                mass[0] * h[0] + incoming * partial_h(500.0, -4000.0)[1],
                incoming,
                0.0,
            )
        };
        let expected_mass = expected_species.iter().sum::<f64>();
        near(result.carrier.projection.field.flux_x(1, 0, 0), rate, 1e-12);
        near(result.carrier.mass_kg[0], expected_mass, 1e-12);
        near(
            result.carrier.specific_enthalpy_j_kg[0],
            total_h / expected_mass,
            1e-6,
        );
        for k in 0..2 {
            near(
                result.species.species_mass_kg[0][k],
                expected_species[k],
                1e-12,
            );
            near(
                result.species.report.outgoing_mass_kg[k],
                outgoing_mass * constituents[0][k] / mass[0],
                1e-12,
            );
        }
        near(
            result.state.temperatures_k()[0],
            temperature(total_h, expected_species, -4000.0),
            1e-8,
        );
        near(result.carrier.report.incoming_mass_kg, incoming_mass, 1e-12);
        near(result.carrier.report.outgoing_mass_kg, outgoing_mass, 1e-12);
        near(
            result.carrier.report.enthalpy_change_j,
            total_h - mass[0] * h[0],
            1e-6,
        );
        near(
            result.carrier.report.carried_enthalpy_j,
            total_h - mass[0] * h[0],
            1e-6,
        );
        assert_eq!(result.carrier.report.mass_routes, 1);
        assert_eq!(result.species.report.source_incoming_mass_kg, [0.0; 2]);
        near(
            result.species.report.boundary_incoming_mass_kg[0],
            0.0,
            1e-12,
        );
        near(
            result.species.report.boundary_incoming_mass_kg[1],
            incoming_mass,
            1e-12,
        );
        assert_eq!(result.report.source_evaluations, 2);
        assert!(result.report.absolute_enthalpy_defect_j <= ENERGY_BUDGET);
    }
    assert_eq!(flow, before);
    assert_eq!(mass, [2.0]);
    assert_eq!(constituents, [[1.5, 0.5]]);
}
