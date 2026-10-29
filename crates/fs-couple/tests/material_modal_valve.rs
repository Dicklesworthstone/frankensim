//! Original anisotropic tetrahedral viscosity through the existing modal
//! aperture. Synthetic geometry/laws isolate coupling, not measured reed data.
use fs_couple::bernoulli_aperture::dynamic::ApertureDrive;
use fs_couple::bernoulli_aperture::modal::{ModalAperture, solid::MaterialSolidOptions};
use fs_couple::bernoulli_aperture::paired::{BladePropagator, BladeState, PairedApertureLimits};
use fs_dcontact::Obstacle;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::gas::{GasSpec, GasState};
use fs_material::visco::tensor::{Mandel6, TensorProny};
use fs_modal::ModalDampingBudget;
use fs_solid::linear3::{
    TetAssemblyBudget,
    visco::{TetViscoField, TetViscoMaterial, ViscoStateBudget},
};
use fs_solid::reduce::{
    OrificeSpec, ValveCardRequest, ValveElementOrder,
    visco::{Controls, MaterialModalValve, reduce_modal_pressure_loaded_valve},
};

const NODES: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const FACE: [usize; 3] = [1, 2, 3];
const PRESSURE: [[f64; 3]; 4] = [[0.; 3], [0.; 3], [0., 0.002, 0.], [0.; 3]];
const FORCE: [[f64; 3]; 4] = [[0.; 3], [0.; 3], [0., 1., 0.], [0., 0., 0.5]];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 83,
                kernel_id: 83,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn diagonal(values: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { values[i] } else { 0.0 }))
}
fn image(pressure: bool, coupled: bool) -> MaterialModalValve {
    image_with_strain_limit(pressure, coupled, None)
}

fn image_with_strain_limit(pressure: bool, coupled: bool, strain_limit: Option<f64>) -> MaterialModalValve {
    with_cx(|cx| {
        let stiffness = diagonal([80.0, 100.1, 100.0, 40.0, 45.0, 50.0]);
        let mut viscosity = diagonal([0.08, 0.1, 0.1, 0.04, 0.045, 0.05]);
        if coupled {
            viscosity[1][2] = 0.09;
            viscosity[2][1] = 0.09;
        }
        let law = if let Some(limit) = strain_limit {
            // An elastic bounded law isolates source-domain admission from
            // the separate full-viscosity assertions in the existing fixture.
            TensorProny::new_mandel(stiffness, vec![], limit).unwrap()
                .with_band((0.0, 10_000.0)).unwrap()
        } else {
            TensorProny::new(stiffness, Some(viscosity), vec![], (0.0, 10_000.0)).unwrap()
        };
        let material = TetViscoMaterial {
            law,
            density_kg_m3: 1.0,
            material_state_identity: [83; 32],
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
            window_hz: (0.1, 100.0),
            retained_compliance_floor: 0.8,
            loss_factor: 0.0,
            source_id: "synthetic/complete-tetrahedral-modal-aperture",
        };
        let controls = Controls {
            reference_omega_rad_s: 10.0,
            maximum_modes: 8,
            assembly: TetAssemblyBudget::standard(),
            internal_variables: ViscoStateBudget {
                max_internal_variables: 1000,
            },
            maximum_storage_sample_error: 0.1,
            maximum_loss_sample_error: 0.02,
            maximum_modal_loss_factor: 0.2,
            damping: ModalDampingBudget {
                max_rel_decay_error: 0.001,
                max_freq_error_cents: 0.01,
                max_cluster_size: 8,
            },
        };
        reduce_modal_pressure_loaded_valve(
            &request,
            if pressure { &PRESSURE } else { &[[0.; 3]; 4] },
            &controls,
            cx,
        )
        .unwrap()
    })
}
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn options(dt: f64) -> MaterialSolidOptions {
    MaterialSolidOptions {
        impedance_pa_s_m3: 1e4,
        time_step_s: dt,
        max_steps: 1024,
        limits: PairedApertureLimits {
            maximum_abs_displacement_m: 0.01,
            maximum_abs_velocity_m_s: 1.0,
            maximum_energy_j: 1.0,
            max_evaluations: 256,
            velocity_tolerance_m_s: 1e-13,
            relative_tolerance: 1e-11,
        },
    }
}
fn contact() -> Obstacle {
    Obstacle::new(
        vec![-1.],
        1,
        1,
        vec![0.],
        vec![1.],
        1e8,
        2.,
        "synthetic fixed mean-gap obstacle".into(),
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-25),
        "{actual:e} != {expected:e}"
    );
}
fn zero_drive() -> ApertureDrive {
    ApertureDrive {
        upstream_pressure_pa: 0.0,
        incoming_pressure_pa: 0.0,
        body_flow_m3_s: 0.0,
    }
}

#[test]
fn material_adapter_keeps_dark_ports_complete_c_and_source_virtual_work() {
    let source = image(true, true);
    let source_identity = source.identity();
    let mass = source.total_mass_kg();
    let root = mass.sqrt();
    let c = source.approximation().modal_viscosity_per_s.clone();
    let state = vec![
        BladeState {
            displacement_m: 2e-6,
            velocity_m_s: 3e-4,
        },
        BladeState {
            displacement_m: -1e-6,
            velocity_m_s: -2e-4,
        },
    ];
    let gate = CancelGate::new_clock_free();
    let aperture = ModalAperture::from_material_solid(
        source,
        [Some(&FORCE), None],
        gas(),
        options(0.001),
        state.clone(),
        contact(),
        &gate,
    )
    .unwrap();
    assert_eq!(aperture.blade_propagator(), BladePropagator::ExactZoh);
    let solid = aperture.material_solid().unwrap();
    assert_eq!(solid.source().identity(), source_identity);
    assert_eq!(aperture.spec().modes.len(), 2);
    let dark = aperture
        .spec()
        .modes
        .iter()
        .position(|m| m.gap_projection == 0.0)
        .unwrap();
    assert_ne!(aperture.spec().modes[dark].pressure_area_m2, 0.0);
    assert_ne!(aperture.spec().modes[dark].force_projection[0], 0.0);
    for (i, mode) in aperture.spec().modes.iter().enumerate() {
        assert_eq!(mode.mass_kg, mass);
        assert_eq!(mode.damping_n_s_m, 0.0);
        close(
            mode.stiffness_n_m,
            mass * solid.source().modes()[i].lambda,
            1e-14,
        );
        close(
            mode.pressure_area_m2,
            root * solid.source().pressure_projection_m2_per_sqrt_kg()[i],
            1e-14,
        );
    }
    for (actual, expected) in aperture.viscosity_matrix_n_s_m().unwrap().iter().zip(c) {
        close(*actual, mass * expected, 1e-14);
    }
    let (nodal_displacement, nodal_velocity) = aperture.solid_nodal_motion(&gate).unwrap();
    for (nodal, coordinate) in [(&nodal_displacement, 0), (&nodal_velocity, 1)] {
        let modal: Vec<_> = state
            .iter()
            .map(|s| {
                if coordinate == 0 {
                    s.displacement_m
                } else {
                    s.velocity_m_s
                }
            })
            .collect();
        let physical_force_work: f64 = nodal
            .iter()
            .flatten()
            .zip(FORCE.iter().flatten())
            .map(|(v, f)| v * f)
            .sum();
        let modal_force_work: f64 = modal
            .iter()
            .zip(&aperture.spec().modes)
            .map(|(v, m)| v * m.force_projection[0])
            .sum();
        close(modal_force_work, physical_force_work, 1e-13);
        let physical_pressure_work: f64 = nodal
            .iter()
            .flatten()
            .zip(PRESSURE.iter().flatten())
            .map(|(v, p)| v * p)
            .sum();
        let modal_pressure_work: f64 = modal
            .iter()
            .zip(&aperture.spec().modes)
            .map(|(v, m)| v * m.pressure_area_m2)
            .sum();
        close(modal_pressure_work, physical_pressure_work, 1e-13);
    }
}

#[test]
fn actual_full_material_loss_excites_dark_mode_and_closes_free_decay_energy() {
    let source = image(false, true);
    let dark = source
        .opening_projection_per_sqrt_kg()
        .iter()
        .position(|g| *g == 0.0)
        .unwrap();
    let bright = 1 - dark;
    let mut state = vec![BladeState::default(); 2];
    state[bright].velocity_m_s = 1e-3;
    let gate = CancelGate::new_clock_free();
    let build = |source, state, dt| {
        ModalAperture::from_material_solid(
            source,
            [None, None],
            gas(),
            options(dt),
            state,
            contact(),
            &gate,
        )
        .unwrap()
    };
    let mut full = build(source, state.clone(), 0.001);
    let mut uncoupled = build(image(false, false), state.clone(), 0.001);
    let mut refined = build(image(false, true), state, 0.0005);
    let initial = full.stored_energy_j();
    let mut heat = 0.0;
    for _ in 0..200 {
        let frame = full.step(zero_drive(), [0.0; 2]).unwrap();
        heat += frame.aperture.dissipated_energy_j;
        assert!(frame.aperture.dissipated_energy_j >= 0.0);
        assert!(frame.aperture.balance_residual_j().abs() < 1e-9 * initial);
        uncoupled.step(zero_drive(), [0.0; 2]).unwrap();
        refined.step(zero_drive(), [0.0; 2]).unwrap();
        refined.step(zero_drive(), [0.0; 2]).unwrap();
    }
    assert!(full.modal_states()[dark].velocity_m_s.abs() > 1e-6);
    assert_eq!(uncoupled.modal_states()[dark], BladeState::default());
    assert!(heat > 0.01 * initial);
    close(full.stored_energy_j() + heat, initial, 1e-8);
    for (coarse, fine) in full.modal_states().iter().zip(refined.modal_states()) {
        close(coarse.displacement_m, fine.displacement_m, 1e-8);
        close(coarse.velocity_m_s, fine.velocity_m_s, 1e-8);
    }
}

#[test]
fn material_adapter_requires_all_states_and_preserves_cancellation() {
    let gate = CancelGate::new_clock_free();
    assert!(
        ModalAperture::from_material_solid(
            image(false, true),
            [None, None],
            gas(),
            options(0.001),
            vec![BladeState::default()],
            contact(),
            &gate,
        )
        .is_err()
    );
    let source = image(false, true);
    gate.request();
    assert!(matches!(
        ModalAperture::from_material_solid(
            source,
            [None, None],
            gas(),
            options(0.001),
            vec![BladeState::default(); 2],
            contact(),
            &gate,
        ),
        Err(fs_couple::acoustic_realize::AcousticRealizeError::Cancelled)
    ));
}

#[path = "material_modal_valve/hereditary.rs"]
mod hereditary;

#[test]
fn original_solid_strain_limits_gate_cold_and_actual_forced_states_atomically() {
    let gate = CancelGate::new_clock_free();
    let build = |state| ModalAperture::from_material_solid(
        image_with_strain_limit(false, false, Some(1e-10)),
        [Some(&FORCE), None], gas(), options(0.001), state, contact(), &gate,
    );
    // These modal amplitudes are inside the independent 1 cm displacement
    // allowance, but violate the actual element source's much smaller strain.
    assert!(build(vec![BladeState { displacement_m: 1e-6, velocity_m_s: 0.0 }; 2]).is_err());
    let mut actual = build(vec![BladeState::default(); 2]).unwrap();
    let mut control = build(vec![BladeState::default(); 2]).unwrap();
    let state = actual.modal_states().to_vec();
    let energy = actual.stored_energy_j();
    let error = actual.step(zero_drive(), [0.01, 0.0]).unwrap_err();
    assert!(error.to_string().contains("small-strain"), "{error}");
    assert_eq!(actual.modal_states(), state.as_slice());
    assert_eq!(actual.stored_energy_j().to_bits(), energy.to_bits());
    actual.step(zero_drive(), [0.0; 2]).unwrap();
    control.step(zero_drive(), [0.0; 2]).unwrap();
    assert_eq!(actual.modal_states(), control.modal_states());
    assert_eq!(actual.stored_energy_j().to_bits(), control.stored_energy_j().to_bits());

    let image = image_with_strain_limit(false, false, Some(1e-10));
    let mut displacement = vec![0.0; 12];
    displacement[7] = 1e-6; // yy strain equals this value in the unit tetrahedron.
    assert!(image.source().validate_displacement(&displacement, || false).is_err());
    displacement[7] = 5e-11;
    image.source().validate_displacement(&displacement, || false).unwrap();
    assert!(image.validate_modal_displacement(&[0.0; 2], || true).is_err());
    image.validate_modal_displacement(&[0.0; 2], || false).unwrap();
}
