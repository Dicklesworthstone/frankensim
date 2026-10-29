//! Synthetic source laws isolate material-to-valve dynamics, not measured reeds.
use fs_exec::Cx;
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_modal::ModalDampingBudget;
use fs_solid::linear3::{
    TetAssemblyBudget,
    visco::{TetViscoField, TetViscoMaterial, ViscoStateBudget},
};
use fs_solid::reduce::{
    OrificeSpec, SingleModeValveAdmission, ValveCard, ValveCardRequest, ValveElementOrder,
    visco::{Controls, reduce_pressure_loaded_valve},
};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const FIXED: [usize; 11] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
const FACE: [usize; 3] = [1, 2, 3];
const PRESSURE: [[f64; 3]; 4] = [[0.0; 3], [0.0; 3], [0.0; 3], [0.0, 0.0, 0.002]];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 80,
                kernel_id: 80,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}
fn diagonal(values: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { values[i] } else { 0.0 }))
}
fn material(beta_s: f64, branch_fraction: f64, tau_s: f64) -> TetViscoMaterial {
    let stiffness = diagonal([80.0, 100.1, 100.0, 40.0, 45.0, 50.0]);
    let branches = if branch_fraction == 0.0 {
        vec![]
    } else {
        vec![
            PronyBranch::new(stiffness.map(|row| row.map(|x| x * branch_fraction)), tau_s).unwrap(),
        ]
    };
    TetViscoMaterial {
        law: TensorProny::new(
            stiffness,
            Some(stiffness.map(|row| row.map(|x| x * beta_s))),
            branches,
            (0.0, 10_000.0),
        )
        .unwrap(),
        density_kg_m3: 1.0,
        material_state_identity: [80; 32],
    }
}
fn controls() -> Controls {
    Controls {
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
    }
}
fn request(field: TetViscoField<'_>) -> ValveCardRequest<'_, TetViscoField<'_>> {
    ValveCardRequest {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        element_order: ValveElementOrder::Tet4,
        fixed_dofs: &FIXED,
        material: field,
        orifice: OrificeSpec {
            face_nodes: &FACE,
            opening_axis: [0.0, 0.0, 1.0],
            width_axis: [1.0, 0.0, 0.0],
            opposing_plane_offset_m: -0.001,
        },
        window_hz: (0.1, 100.0),
        retained_compliance_floor: 0.8,
        loss_factor: 0.0,
        source_id: "synthetic/single-axial-tet-valve",
    }
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-25),
        "{actual:e} != {expected:e}"
    );
}
fn admission() -> SingleModeValveAdmission {
    SingleModeValveAdmission {
        minimum_opening_compliance_fraction: 0.8,
        minimum_pressure_compliance_fraction: 0.8,
        maximum_static_opening_relative_error: 0.1,
    }
}
fn decay(card: &ValveCard) -> (f64, f64) {
    let mode = &card.modes[0];
    let system =
        fs_phs::mass_spring_damper(mode.mass_kg, mode.stiffness_n_m, mode.damping_n_s_m).unwrap();
    let mut state = vec![1e-4, 0.0];
    let initial = system.hamiltonian(&state);
    let mut dissipated = 0.0;
    for _ in 0..2000 {
        let step = fs_phs::step(&system, &state, &[0.0], 0.001).unwrap();
        dissipated += step.dissipated;
        state = step.x;
    }
    let final_energy = system.hamiltonian(&state);
    close(final_energy + dissipated, initial, 1e-7);
    (final_energy, dissipated)
}

#[test]
fn source_kelvin_viscosity_changes_actual_decay_and_playing_reed() {
    with_cx(|cx| {
        let low = material(0.001, 0.0, 1.0);
        let high = material(0.005, 0.0, 1.0);
        let a = reduce_pressure_loaded_valve(
            &request(TetViscoField::Uniform(&low)),
            &PRESSURE,
            &controls(),
            cx,
        )
        .unwrap();
        let b = reduce_pressure_loaded_valve(
            &request(TetViscoField::Uniform(&high)),
            &PRESSURE,
            &controls(),
            cx,
        )
        .unwrap();
        let am = &a.valve().card().modes[0];
        let bm = &b.valve().card().modes[0];
        // V=1/6, M_zz=rho*V/10 and unit mean opening makes u_tip=3.
        close(am.mass_kg, 0.15, 1e-12);
        close(am.stiffness_n_m, 150.0, 1e-12);
        close(am.damping_n_s_m, 0.15, 1e-12);
        close(bm.mass_kg, am.mass_kg, 1e-12);
        close(bm.stiffness_n_m, am.stiffness_n_m, 1e-12);
        close(bm.damping_n_s_m, 5.0 * am.damping_n_s_m, 1e-12);
        close(a.valve().effective_areas_m2()[0], 0.006, 1e-12);
        assert_eq!(a.valve().card().loss_factor, 0.0);
        let ar = a.beating_reed_checked(1.0, 0.0, admission()).unwrap();
        let br = b.beating_reed_checked(1.0, 0.0, admission()).unwrap();
        close(
            ar.damping_ratio,
            am.damping_n_s_m / (2.0 * (am.mass_kg * am.stiffness_n_m).sqrt()),
            1e-13,
        );
        close(br.damping_ratio, 5.0 * ar.damping_ratio, 1e-12);
        let (low_energy, low_heat) = decay(a.valve().card());
        let (high_energy, high_heat) = decay(b.valve().card());
        assert!(high_energy < 0.05 * low_energy);
        assert!(high_heat > low_heat);
        assert_ne!(a.valve().identity(), b.valve().identity());
        let decoded =
            ValveCard::from_canonical_bytes(&b.valve().card().to_canonical_bytes()).unwrap();
        let reed = decoded.beating_reed(0.006, 1.0, 0.0).unwrap();
        close(reed.damping_ratio, br.damping_ratio, 1e-14);
    })
}

#[test]
fn prony_data_changes_frequency_and_loss_without_relabeling_equilibrium_compliance() {
    with_cx(|cx| {
        let first = material(0.001, 0.05, 0.02);
        let second = material(0.001, 0.05, 0.03);
        let options = controls();
        let a = reduce_pressure_loaded_valve(
            &request(TetViscoField::Uniform(&first)),
            &PRESSURE,
            &options,
            cx,
        )
        .unwrap();
        let b = reduce_pressure_loaded_valve(
            &request(TetViscoField::Uniform(&second)),
            &PRESSURE,
            &options,
            cx,
        )
        .unwrap();
        for (law, valve) in [(&first, &a), (&second, &b)] {
            let (storage, loss) = law
                .law
                .complex_modulus(options.reference_omega_rad_s)
                .unwrap();
            let mode = &valve.valve().card().modes[0];
            close(mode.stiffness_n_m, 1.5 * storage[2][2], 1e-12);
            close(
                mode.damping_n_s_m,
                1.5 * loss[2][2] / options.reference_omega_rad_s,
                1e-12,
            );
            close(
                mode.frequency_hz,
                (10.0 * storage[2][2]).sqrt() / std::f64::consts::TAU,
                1e-12,
            );
            // All degrees are retained, but source storage at omega_ref is
            // stiffer than equilibrium. The static comparator MUST still see it.
            close(
                valve.valve().card().retained_compliance_fraction,
                100.0 / storage[2][2],
                1e-12,
            );
            close(
                valve
                    .valve()
                    .pressure_compliance()
                    .exact_mean_force_m_per_n(),
                1.0 / 150.0,
                1e-12,
            );
            assert!(valve.approximation().maximum_storage_sample_error > 1e-5);
            assert!(valve.approximation().maximum_loss_sample_error > 1e-5);
        }
        assert!(b.valve().card().modes[0].frequency_hz > a.valve().card().modes[0].frequency_hz);
        assert!(b.valve().card().modes[0].damping_n_s_m > a.valve().card().modes[0].damping_n_s_m);
        assert_ne!(
            a.valve().card().material_identity,
            b.valve().card().material_identity
        );
        let expanded = [first.clone()];
        let explicit = reduce_pressure_loaded_valve(
            &request(TetViscoField::PerElement(&expanded)),
            &PRESSURE,
            &options,
            cx,
        )
        .unwrap();
        assert_eq!(a.valve().card(), explicit.valve().card());
        assert_eq!(a.valve().identity(), explicit.valve().identity());
    })
}

#[test]
fn source_bands_double_count_and_dispersion_policy_are_enforced() {
    with_cx(|cx| {
        let m = material(0.001, 0.05, 0.02);
        let mut req = request(TetViscoField::Uniform(&m));
        req.loss_factor = 0.01;
        assert!(reduce_pressure_loaded_valve(&req, &PRESSURE, &controls(), cx).is_err());
        req.loss_factor = 0.0;
        let mut strict = controls();
        strict.maximum_storage_sample_error = 1e-5;
        assert!(reduce_pressure_loaded_valve(&req, &PRESSURE, &strict, cx).is_err());
        strict = controls();
        strict.maximum_loss_sample_error = 1e-5;
        assert!(reduce_pressure_loaded_valve(&req, &PRESSURE, &strict, cx).is_err());
        for band in [(1.0, 1000.0), (0.0, 5.0), (0.0, 20.0)] {
            let mut limited = m.clone();
            limited.law = limited.law.with_band(band).unwrap();
            assert!(
                reduce_pressure_loaded_valve(
                    &request(TetViscoField::Uniform(&limited)),
                    &PRESSURE,
                    &controls(),
                    cx
                )
                .is_err()
            );
        }
    })
}

#[test]
fn nonproportional_near_degenerate_material_loss_cannot_be_exported_as_independent_modes() {
    with_cx(|cx| {
        let mut m = material(0.001, 0.0, 1.0);
        let mut viscosity = diagonal([0.08, 0.1, 0.1, 0.04, 0.045, 0.05]);
        viscosity[1][2] = 0.09;
        viscosity[2][1] = 0.09;
        m.law = TensorProny::new(*m.law.c_eq(), Some(viscosity), vec![], (0.0, 10_000.0)).unwrap();
        let fixed: Vec<usize> = (0..12).filter(|&d| d != 7 && d != 11).collect();
        let mut req = request(TetViscoField::Uniform(&m));
        req.fixed_dofs = &fixed;
        req.orifice.opening_axis = [0.0, 1.0, 1.0];
        let pressure = [[0.0; 3], [0.0; 3], [0.0, 0.002, 0.0], [0.0, 0.0, 0.002]];
        assert!(reduce_pressure_loaded_valve(&req, &pressure, &controls(), cx).is_err());
        // Same elastic material and geometry, proportional viscosity: admitted.
        let uncoupled = material(0.001, 0.0, 1.0);
        req.material = TetViscoField::Uniform(&uncoupled);
        assert!(reduce_pressure_loaded_valve(&req, &pressure, &controls(), cx).is_ok());
    })
}


#[test]
fn complete_image_retains_dark_pressure_mode_and_original_material_cross_loss() {
    use fs_modal::ModalDampingBlock;
    use fs_solid::reduce::visco::reduce_modal_pressure_loaded_valve;
    with_cx(|cx| {
        let mut m = material(0.001, 0.0, 1.0);
        let mut viscosity = diagonal([0.08, 0.1, 0.1, 0.04, 0.045, 0.05]);
        viscosity[1][2] = 0.09;
        viscosity[2][1] = 0.09;
        m.law = TensorProny::new(*m.law.c_eq(), Some(viscosity), vec![], (0.0, 10_000.0)).unwrap();
        let fixed: Vec<_> = (0..12).filter(|&d| d != 7 && d != 11).collect();
        let mut req = request(TetViscoField::Uniform(&m));
        req.fixed_dofs = &fixed;
        // Mean opening observes node 3's z coordinate. Pressure acts solely
        // on node 2's y coordinate, which the compatibility card would omit.
        let pressure = [[0.0; 3], [0.0; 3], [0.0, 0.002, 0.0], [0.0; 3]];
        assert!(reduce_pressure_loaded_valve(&req, &pressure, &controls(), cx).is_err());
        let image = reduce_modal_pressure_loaded_valve(&req, &pressure, &controls(), cx).unwrap();
        assert_eq!(image.modes().len(), 2);
        assert!(image.approximation().damping.blocks.iter()
            .any(|block| matches!(block, ModalDampingBlock::Full { .. })));
        let dark = image.opening_projection_per_sqrt_kg().iter()
            .position(|g| *g == 0.0).expect("exact zero-opening mode retained");
        assert_ne!(image.pressure_projection_m2_per_sqrt_kg()[dark], 0.0);
        close(image.compliance().retained_volume_m3_per_pa,
            image.compliance().exact_volume_m3_per_pa, 1e-12);
        close(image.compliance().retained_opening_force_m_per_n,
            image.compliance().exact_opening_force_m_per_n, 1e-12);
        assert_eq!(image.compliance().retained_opening_m_per_pa, 0.0);
        let c = &image.approximation().modal_viscosity_per_s;
        close(c[1].abs(), 0.9, 1e-12);
        assert_eq!(c[1], c[2]);

        let gate = fs_exec::CancelGate::new_clock_free();
        let q = [0.03, -0.04];
        let nodal = image.expand_nodal(&q, &gate).unwrap();
        let force = [[0.0; 3], [0.0; 3], [0.0, 0.7, 0.0], [0.0, 0.0, -0.2]];
        let projected = image.project_nodal_force(&force, &gate).unwrap();
        let nodal_work: f64 = nodal.iter().flatten().zip(force.iter().flatten())
            .map(|(x, f)| x * f).sum();
        let modal_work: f64 = q.iter().zip(&projected).map(|(x, f)| x * f).sum();
        close(modal_work, nodal_work, 1e-13);
        assert_eq!(image.project_nodal_force(&pressure, &gate).unwrap(),
            image.pressure_projection_m2_per_sqrt_kg());

        // Compare to the actual element owner's assembled source viscosity,
        // not a second copy of the modal projection formula.
        let reference = image.source().complex_stiffness(controls().reference_omega_rad_s, cx).unwrap();
        let free_velocity: Vec<_> = image.source().equilibrium().free_dofs.iter()
            .map(|&d| nodal[d / 3][d % 3]).collect();
        let mut source_force = vec![0.0; free_velocity.len()];
        reference.imag.spmv(&free_velocity, &mut source_force);
        let source_power: f64 = free_velocity.iter().zip(&source_force)
            .map(|(v, f)| v * f / reference.omega_rad_s).sum();
        let modal_power: f64 = (0..2).map(|i| (0..2)
            .map(|j| q[i] * c[2*i+j] * q[j]).sum::<f64>()).sum();
        close(modal_power, source_power, 1e-12);
        assert!(modal_power > 0.0);
        assert!(image.project_nodal_force(&force[..3], &gate).is_err());
        assert!(image.expand_nodal(&[0.0], &gate).is_err());
        gate.request();
        assert!(image.expand_nodal(&q, &gate).is_err());
    });
}

#[test]
fn complete_material_image_keeps_source_admission_and_legacy_identity() {
    use fs_solid::reduce::visco::reduce_modal_pressure_loaded_valve;
    with_cx(|cx| {
        let m = material(0.001, 0.05, 0.02);
        let mut req = request(TetViscoField::Uniform(&m));
        let image = reduce_modal_pressure_loaded_valve(&req, &PRESSURE, &controls(), cx).unwrap();
        let legacy = reduce_pressure_loaded_valve(&req, &PRESSURE, &controls(), cx).unwrap();
        assert_eq!(image.modes().len(), legacy.valve().source_modes().len());
        for (complete, compatible) in image.modes().iter().zip(legacy.valve().source_modes()) {
            assert_eq!(complete.lambda, compatible.lambda);
            assert_eq!(complete.phi, compatible.phi);
        }
        assert_eq!(image.approximation(), legacy.approximation());
        let first_identity = image.identity();
        assert_eq!(first_identity,
            reduce_modal_pressure_loaded_valve(&req, &PRESSURE, &controls(), cx).unwrap().identity());
        req.loss_factor = 0.1;
        assert!(reduce_modal_pressure_loaded_valve(&req, &PRESSURE, &controls(), cx).is_err());
        req.loss_factor = 0.0;
        let mut strict = controls();
        strict.maximum_storage_sample_error = 1e-5;
        assert!(reduce_modal_pressure_loaded_valve(&req, &PRESSURE, &strict, cx).is_err());
        strict = controls();
        strict.maximum_loss_sample_error = 1e-5;
        assert!(reduce_modal_pressure_loaded_valve(&req, &PRESSURE, &strict, cx).is_err());
        let changed = [[0.0; 3], [0.0; 3], [0.0; 3], [0.0, 0.0, 0.004]];
        assert_ne!(first_identity,
            reduce_modal_pressure_loaded_valve(&req, &changed, &controls(), cx).unwrap().identity());
    });
}

#[path = "visco_valve/hereditary.rs"]
mod hereditary;
