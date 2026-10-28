//! Full 3D tensor material reaches the SAME original reed trajectory and ports.
use super::*;
use fs_couple::bernoulli_aperture::modal::relaxation::{
    InitialModalMemory, TensorModalRelaxationLimits,
};
use fs_couple::bernoulli_aperture::paired::BladePropagator;
use fs_material::visco_tensor::{TensorProny, TensorPronyBranch};
use fs_material::{IsotropicElastic, Tangent6};

fn bare() -> ModalAperture {
    // A deliberately coarse two-triangle manufactured specimen. Its original
    // curvature/eigenbasis is retained; this is not measured cane geometry.
    let mesh = PlateMesh::rectangle(0.025, 0.01, 1, 1);
    let supports = mesh
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| p.0 == 0.0)
        .map(|(i, _)| i)
        .collect();
    let slit: Vec<_> = mesh
        .boundary_edges()
        .into_iter()
        .filter(|&(a, b)| mesh.nodes[a].0 == 0.025 && mesh.nodes[b].0 == 0.025)
        .map(|(a, b)| [a, b])
        .collect();
    let mut force = vec![0.0; mesh.nodes.len()];
    force[slit[0][0]] = 1.0;
    let elements = (0..mesh.tris.len()).collect();
    let chart = compile_plate_material_chart(
        mesh,
        supports,
        vec![PlateRegionMaterial {
            region: PlateRegion {
                name: "original tensor blade".into(),
                triangle_indices: elements,
            },
            material: source(900.0, 0.0),
            model: PlateMaterialModel::Isotropic,
            thickness_constraint: PlateThicknessConstraint::FixedThickness(0.0003),
            material_angle_rad: 0.0,
        }],
    )
    .unwrap();
    let options = ModalPlateOptions {
        reduction: PlateChartRadiation {
            assembly: AssemblyOptions {
                pretension: 0.0,
                support: EdgeSupport::Clamped,
            },
            eigenvalue_window: (1.0, (std::f64::consts::TAU * 15_000.0).powi(2)),
            n_modes: 2,
            damping_ratio: 0.0,
            unit_force_weights: force,
        },
        slit_edges: slit,
        maximum_trace_slope: 0.1,
        maximum_slit_deviation_m: 0.001,
        maximum_angular_step: 1.0,
        max_nodes: 16,
        max_triangles: 8,
        max_projection_terms: 1_000_000,
    };
    ModalAperture::from_material_plates_with_propagator(
        vec![(chart, options)],
        0.0004,
        super::super::network::gas(),
        1e6,
        1.0 / 192_000.0,
        1024,
        limits(),
        vec![BladeState::default(); 2],
        contact(),
        BladePropagator::ExactZoh,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn law(strength: f64, strain_limit: f64) -> TensorProny {
    let mut factor: Tangent6 = [[0.0; 6]; 6];
    factor[0] = [10000.0 * strength, -2000.0, 12000.0, 4000.0, 1000.0, -700.0];
    TensorProny::from_tensor_shear(
        IsotropicElastic::new(4e9, 0.3, 0.1).unwrap().stiffness(),
        vec![TensorPronyBranch::from_gram_rows(factor, 0.001).unwrap()],
        strain_limit,
    )
    .unwrap()
    .with_band((0.0, 100000.0))
    .unwrap()
}
fn capacity() -> TensorModalRelaxationLimits {
    TensorModalRelaxationLimits {
        max_material_points: 6,
        max_memory_states: 36,
        max_projection_terms: 16_000_000,
    }
}
fn build(strength: f64, strain_limit: f64) -> ModalAperture {
    bare()
        .with_tensor_material_relaxation(
            vec![vec![law(strength, strain_limit)]],
            InitialModalMemory::Relaxed,
            capacity(),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
}

#[test]
fn full_tensor_memory_changes_material_reed_pressure_with_original_force_and_loss_audits() {
    let mut a = build(1.0, 0.1);
    let mut b = build(3.0, 0.1);
    let original = a.material_plates()[0].source().specimen_identity();
    assert_eq!(
        original,
        b.material_plates()[0].source().specimen_identity()
    );
    assert_eq!(
        a.spec().modes,
        b.spec().modes,
        "only the full tensor memory changed"
    );
    let material = a.material_relaxation().unwrap();
    assert_eq!(material.tensor_material_points().len(), 6);
    assert_eq!(
        material.memory_sqrt_j()[0].len(),
        36,
        "six original rows per point survive the rank-one law"
    );
    assert!(
        material.material_spectra().is_empty(),
        "a full tensor must not advertise a scalar Young spectrum"
    );
    let mut changed = false;
    let mut dissipated = 0.0;
    for sample in 0..512 {
        let input = drive(if sample < 128 { 30.0 } else { 0.0 });
        let x = a.step(input, [0.0; 2]).unwrap();
        let y = b.step(input, [0.0; 2]).unwrap();
        changed |= (x.aperture.outgoing_pressure_pa - y.aperture.outgoing_pressure_pa).abs() > 1e-8;
        assert!(x.work_certificate.residual <= x.work_certificate.effective_budget);
        assert!(x.flow_certificate.residual <= x.flow_certificate.effective_budget);
        assert!(x.momentum_residual_n.iter().all(|v| v.abs() < 1e-7));
        dissipated += x.aperture.relaxation_loss_j;
    }
    assert!(changed);
    assert!(dissipated > 0.0);
    let material = a.material_relaxation().unwrap();
    let z = &material.memory_sqrt_j()[0];
    for point in material.tensor_material_points() {
        let q: Vec<_> = a.modal_states()[point.mode_range()]
            .iter()
            .map(|s| s.displacement_m)
            .collect();
        let normal = point
            .image()
            .normal_stress_pa(&q, &z[point.memory_range()])
            .unwrap();
        assert!(normal.iter().all(|s| s.abs() < 1e-5), "{normal:?}");
        assert!(point.surface_strain(a.modal_states(), z).is_ok());
    }
}

#[test]
fn tensor_reed_refuses_missing_equilibrium_and_history_capacity_without_dropping_modes() {
    let gate = CancelGate::new_clock_free();
    let mut cap = capacity();
    cap.max_memory_states = 35;
    assert!(
        bare()
            .with_tensor_material_relaxation(
                vec![vec![law(1.0, 0.1)]],
                InitialModalMemory::Relaxed,
                cap,
                &gate
            )
            .is_err()
    );
    cap = capacity();
    cap.max_material_points = 5;
    assert!(
        bare()
            .with_tensor_material_relaxation(
                vec![vec![law(1.0, 0.1)]],
                InitialModalMemory::Relaxed,
                cap,
                &gate
            )
            .is_err()
    );
    cap = capacity();
    cap.max_projection_terms = 1;
    assert!(
        bare()
            .with_tensor_material_relaxation(
                vec![vec![law(1.0, 0.1)]],
                InitialModalMemory::Relaxed,
                cap,
                &gate
            )
            .is_err()
    );
    let mismatch = TensorProny::from_tensor_shear(
        IsotropicElastic::new(5e9, 0.3, 0.1).unwrap().stiffness(),
        Vec::new(),
        0.1,
    )
    .unwrap();
    assert!(
        bare()
            .with_tensor_material_relaxation(
                vec![vec![mismatch]],
                InitialModalMemory::Relaxed,
                capacity(),
                &gate
            )
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        bare()
            .with_tensor_material_relaxation(
                vec![vec![law(1.0, 0.1)]],
                InitialModalMemory::Relaxed,
                capacity(),
                &cancelled
            )
            .is_err()
    );
}

#[test]
fn physical_surface_strain_refusal_preserves_every_original_reed_history() {
    let mut a = build(1.0, 1e-16);
    let state = a.modal_states().to_vec();
    let memory = a.material_relaxation().unwrap().memory_sqrt_j().to_vec();
    let energy = a.stored_energy_j();
    assert!(a.step(drive(30.0), [0.0; 2]).is_err());
    assert_eq!(a.accepted_steps(), 0);
    assert_eq!(a.modal_states(), state);
    assert_eq!(a.material_relaxation().unwrap().memory_sqrt_j(), memory);
    assert_eq!(a.stored_energy_j(), energy);
    a.step(drive(0.0), [0.0; 2]).unwrap();
}
