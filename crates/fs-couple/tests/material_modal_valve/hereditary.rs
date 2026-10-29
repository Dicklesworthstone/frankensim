use super::*;
use fs_couple::bernoulli_aperture::modal::{
    relaxation::InitialModalMemory, solid::HereditarySolidOptions,
};
use fs_material::visco::tensor::PronyBranch;
use fs_solid::linear3::visco::SolidProjectionBudget;
use fs_solid::reduce::visco::{
    HereditaryModalValve, HereditaryValveControls, reduce_hereditary_pressure_loaded_valve,
};

fn image(pressure: bool, strain_limit: Option<f64>) -> HereditaryModalValve {
    with_cx(|cx| {
        let stiffness = diagonal([80., 100.1, 100., 40., 45., 50.]);
        let mut viscosity = stiffness.map(|r| r.map(|x| 0.001 * x));
        viscosity[1][2] = 0.09;
        viscosity[2][1] = 0.09;
        let mut factor = [[0.; 6]; 6];
        factor[0] = [0., 7., -9., 0., 0., 0.];
        factor[1] = [0., 0.2, 0.1, 0., 0., 0.];
        let branches = vec![PronyBranch::from_gram_rows(factor, 0.02).unwrap()];
        let law = match strain_limit {
            Some(limit) => TensorProny::new_mandel(stiffness, branches, limit).unwrap(),
            None => TensorProny::new(stiffness, Some(viscosity), branches, (0., 10_000.)).unwrap(),
        };
        let material = TetViscoMaterial {
            law,
            density_kg_m3: 1.,
            material_state_identity: [95; 32],
        };
        let fixed: Vec<_> = (0..12).filter(|&d| d != 7 && d != 11).collect();
        let request = ValveCardRequest {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            element_order: ValveElementOrder::Tet4,
            fixed_dofs: &fixed,
            material: TetViscoField::Uniform(&material),
            orifice: OrificeSpec {
                face_nodes: &FACE,
                opening_axis: [0., 0., 1.],
                width_axis: [1., 0., 0.],
                opposing_plane_offset_m: -0.001,
            },
            window_hz: (0.1, 100.),
            retained_compliance_floor: 0.8,
            loss_factor: 0.,
            source_id: "synthetic/hereditary-modal-aperture",
        };
        let controls = HereditaryValveControls {
            assembly: TetAssemblyBudget::standard(),
            internal_variables: ViscoStateBudget {
                max_internal_variables: 1000,
            },
            projection: SolidProjectionBudget {
                max_coordinates: 8,
                max_basis_scalars: 1000,
                max_material_points: 100,
                max_scalar_terms: 1_000_000,
            },
            maximum_factor_scalar_terms: 1_000_000,
            maximum_basis_error: 1e-10,
        };
        reduce_hereditary_pressure_loaded_valve(
            &request,
            if pressure { &PRESSURE } else { &[[0.; 3]; 4] },
            &controls,
            cx,
        )
        .unwrap()
    })
}
fn hereditary_options(dt: f64) -> HereditarySolidOptions {
    HereditarySolidOptions {
        aperture: options(dt),
        maximum_angular_step: 0.3,
    }
}
fn aperture(
    source: HereditaryModalValve,
    state: Vec<BladeState>,
    initial: InitialModalMemory,
    dt: f64,
) -> ModalAperture {
    ModalAperture::from_hereditary_material_solid(
        source,
        [Some(&FORCE), None],
        gas(),
        hereditary_options(dt),
        state,
        initial,
        contact(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

#[test]
fn original_memory_initialization_changes_storage_without_losing_dark_work_ports() {
    let source = image(true, None);
    let state = vec![
        BladeState {
            displacement_m: 2e-5,
            velocity_m_s: 3e-4,
        },
        BladeState {
            displacement_m: -1e-5,
            velocity_m_s: -2e-4,
        },
    ];
    let added = source
        .branches()
        .iter()
        .flat_map(|b| b.rows())
        .map(|row| {
            let strain = row
                .iter()
                .zip(&state)
                .map(|(r, s)| r * s.displacement_m)
                .sum::<f64>();
            0.5 * strain * strain
        })
        .sum::<f64>();
    let identity = source.identity();
    let c = source.viscosity_n_s_m().to_vec();
    let unrelaxed = aperture(source, state.clone(), InitialModalMemory::Unrelaxed, 0.001);
    let relaxed = aperture(image(true, None), state, InitialModalMemory::Relaxed, 0.001);
    assert_eq!(unrelaxed.blade_propagator(), BladePropagator::ExactZoh);
    assert_eq!(
        unrelaxed.hereditary_solid().unwrap().source().identity(),
        identity
    );
    assert_eq!(unrelaxed.viscosity_matrix_n_s_m().unwrap(), c);
    close(
        unrelaxed.stored_energy_j() - relaxed.stored_energy_j(),
        added,
        1e-11,
    );
    let memory = unrelaxed.material_relaxation().unwrap();
    assert!(memory.material_spectra().is_empty());
    assert!(memory.memory_sqrt_j()[0].iter().all(|z| *z == 0.));
    assert!(
        relaxed.material_relaxation().unwrap().memory_sqrt_j()[0]
            .iter()
            .any(|z| *z != 0.)
    );
    let dark = unrelaxed
        .spec()
        .modes
        .iter()
        .position(|m| m.gap_projection == 0.)
        .unwrap();
    assert!(unrelaxed.spec().modes[dark].pressure_area_m2.abs() > 0.);
    assert!(unrelaxed.spec().modes[dark].force_projection[0].abs() > 0.);
    let gate = CancelGate::new_clock_free();
    let (_, velocity) = unrelaxed.solid_nodal_motion(&gate).unwrap();
    let nodal_power = FORCE
        .iter()
        .zip(&velocity)
        .map(|(f, v)| f.iter().zip(v).map(|(a, b)| a * b).sum::<f64>())
        .sum::<f64>();
    let modal_power = unrelaxed
        .spec()
        .modes
        .iter()
        .zip(unrelaxed.modal_states())
        .map(|(m, s)| m.force_projection[0] * s.velocity_m_s)
        .sum::<f64>();
    close(modal_power, nodal_power, 1e-12);
}

#[test]
fn actual_hereditary_ringdown_activates_dark_mode_and_accounts_for_material_heat() {
    let source = image(false, None);
    let opening = source
        .opening_projection_per_sqrt_kg()
        .iter()
        .position(|x| *x != 0.)
        .unwrap();
    let dark = 1 - opening;
    let mut state = vec![BladeState::default(); 2];
    state[opening].displacement_m = 2e-5;
    let mut model = aperture(source, state, InitialModalMemory::Unrelaxed, 0.001);
    let initial = model.stored_energy_j();
    let mut heat = 0.;
    let mut dark_motion = 0.0_f64;
    for _ in 0..500 {
        let frame = model.step(zero_drive(), [0.; 2]).unwrap();
        assert!(frame.aperture.dissipated_energy_j >= 0.);
        heat += frame.aperture.dissipated_energy_j;
        dark_motion = dark_motion.max(model.modal_states()[dark].displacement_m.abs());
    }
    assert!(dark_motion > 1e-7);
    assert!(heat > 0.1 * initial);
    close(model.stored_energy_j() + heat, initial, 2e-8);
    assert!(
        model.material_relaxation().unwrap().memory_sqrt_j()[0]
            .iter()
            .any(|z| *z != 0.)
    );
}

#[test]
fn held_force_response_and_both_loss_channels_are_sample_refinement_invariant() {
    let state = vec![
        BladeState {
            displacement_m: 2e-5,
            velocity_m_s: 3e-4,
        },
        BladeState {
            displacement_m: -1e-5,
            velocity_m_s: -2e-4,
        },
    ];
    let mut coarse = aperture(
        image(false, None),
        state.clone(),
        InitialModalMemory::Relaxed,
        0.002,
    );
    let mut fine = aperture(
        image(false, None),
        state,
        InitialModalMemory::Relaxed,
        0.001,
    );
    let mut loss = [0.; 2];
    for _ in 0..40 {
        loss[0] += coarse
            .step(zero_drive(), [0.0004, 0.])
            .unwrap()
            .aperture
            .dissipated_energy_j;
        for _ in 0..2 {
            loss[1] += fine
                .step(zero_drive(), [0.0004, 0.])
                .unwrap()
                .aperture
                .dissipated_energy_j;
        }
    }
    for (a, b) in coarse.modal_states().iter().zip(fine.modal_states()) {
        close(a.displacement_m, b.displacement_m, 2e-8);
        close(a.velocity_m_s, b.velocity_m_s, 2e-8);
    }
    for (a, b) in coarse.material_relaxation().unwrap().memory_sqrt_j()[0]
        .iter()
        .zip(&fine.material_relaxation().unwrap().memory_sqrt_j()[0])
    {
        close(*a, *b, 2e-8);
    }
    close(loss[0], loss[1], 2e-8);
}

#[test]
fn material_clock_and_endpoint_strain_refusals_preserve_accepted_history() {
    let gate = CancelGate::new_clock_free();
    let mut unresolved = hereditary_options(0.001);
    unresolved.maximum_angular_step = 1e-5;
    assert!(
        ModalAperture::from_hereditary_material_solid(
            image(false, None),
            [None, None],
            gas(),
            unresolved,
            vec![BladeState::default(); 2],
            InitialModalMemory::Relaxed,
            contact(),
            &gate
        )
        .is_err()
    );
    let mut model = aperture(
        image(false, Some(1e-7)),
        vec![BladeState::default(); 2],
        InitialModalMemory::Relaxed,
        0.001,
    );
    let before = model.modal_states().to_vec();
    let memory = model
        .material_relaxation()
        .unwrap()
        .memory_sqrt_j()
        .to_vec();
    assert!(model.step(zero_drive(), [1., 0.]).is_err());
    assert_eq!(model.accepted_steps(), 0);
    assert_eq!(model.modal_states(), before);
    assert_eq!(model.material_relaxation().unwrap().memory_sqrt_j(), memory);
    gate.request();
    assert!(matches!(
        ModalAperture::from_hereditary_material_solid(
            image(false, None),
            [None, None],
            gas(),
            hereditary_options(0.001),
            vec![BladeState::default(); 2],
            InitialModalMemory::Relaxed,
            contact(),
            &gate
        ),
        Err(fs_couple::acoustic_realize::AcousticRealizeError::Cancelled)
    ));
}
