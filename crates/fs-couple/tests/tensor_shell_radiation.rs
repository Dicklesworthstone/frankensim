//! Original tensor shell, full reciprocal BEM loading and coherent receiver fields.
//! These use the actual material, shell, harmonic, boundary and observer owners.

use fs_bem::helmholtz::Medium;
use fs_bem::near_field::{FirstOrder, Options};
use fs_bem::reciprocal_radiation::{Budget, CoherentPortFieldLimits};
use fs_couple::tensor_shell_radiation::{
    TensorShellRadiation, TensorShellRadiationControls, TensorShellRadiationError,
    TensorShellRadiationFrequencyControls,
};
use fs_exec::{CancelGate, Cx};
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;
use fs_modal::harmonic::{HarmonicPortControls, HarmonicResponseControls};
use fs_plate::shell::viscosity::ShellViscosityBudget;
use fs_plate::shell::viscosity::harmonic::{
    TensorPronyShell, TensorPronyShellSection, TensorShellHarmonicError,
};
use fs_plate::shell::{CstDktMass, ShellMesh, ShellSupport};

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(&gate, arena, fs_exec::StreamKey {
            seed: 47, kernel_id: 53, tile: 0, iteration: 0,
        }, fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic), &gate);
    });
}
fn mesh() -> ShellMesh {
    ShellMesh::new(
        vec![[0.0,0.0,0.0], [0.4,0.0,0.0], [0.0,0.3,0.0], [0.0,0.0,0.5]],
        vec![[0,2,1], [0,1,3], [0,3,2], [1,2,3]],
    ).unwrap()
}
fn source(mesh: ShellMesh, cx: &Cx<'_>) -> TensorPronyShell {
    // Explicit numerical polymer-like source and homogeneous water-like fluid;
    // this fixture is not a material-record or spherical-shell benchmark claim.
    let c: Mandel6 = core::array::from_fn(|i| core::array::from_fn(|j| {
        1e6 * if i < 3 && j < 3 {
            4.0 + if i == j { 6.0 } else { 0.0 }
        } else if i == j { 6.0 } else { 0.0 }
    }));
    let scale = |s: f64| c.map(|r| r.map(|x| x*s));
    let law = TensorProny::new_mandel(c, vec![
        PronyBranch::new(scale(0.4),0.001).unwrap(),
        PronyBranch::new(scale(0.2),0.02).unwrap(),
    ], 0.02).unwrap().with_kelvin_viscosity_mandel(scale(1e-5)).unwrap()
        .with_band((0.0,150.0)).unwrap();
    let section = TensorPronyShellSection::new(law,0.01,1200.0).unwrap();
    let sections = vec![section;mesh.tris.len()];
    TensorPronyShell::assemble(mesh, sections, &[0,1,2], ShellSupport::Clamped,
        CstDktMass::Consistent, 0.0, ShellViscosityBudget {
            max_elements:8, max_full_dofs:24, max_dofs:24, max_modes:24,
            max_work:20_000_000,
        }, cx).unwrap()
}
fn medium() -> Medium { Medium { density:1000.0, sound_speed:1480.0 } }
fn geometry_controls() -> TensorShellRadiationControls {
    TensorShellRadiationControls {
        acoustic_refinements:1, maximum_panels:64, maximum_coordinates:6,
        maximum_setup_bytes:1<<20, maximum_setup_work:10_000_000,
    }
}
fn boundary_budget() -> Budget {
    Budget { maximum_unknowns:128, maximum_cache_bytes:16<<20,
        maximum_work:2_000_000_000, angular_refinement:0 }
}
fn frequency_controls() -> TensorShellRadiationFrequencyControls {
    TensorShellRadiationFrequencyControls {
        mechanical:HarmonicResponseControls {
            maximum_coordinates:6, maximum_payload_bytes:16<<20,
            maximum_relative_residual:2e-11,
        },
        ports:HarmonicPortControls {
            maximum_ports:6, maximum_projection_terms:1_000_000,
            maximum_payload_bytes:4<<20,
        },
        maximum_boundary_payload_bytes:16<<20,
        field:CoherentPortFieldLimits { maximum_payload_bytes:1<<20, maximum_work:10_000_000 },
        maximum_radiated_power_work:100_000_000, power_relative_tolerance:1e-8,
    }
}
fn force() -> [C64;6] {
    [C64::new(1e-4,-2e-5), C64::new(-3e-5,1e-5), C64::new(2e-5,4e-5),
        C64::new(1e-6,-2e-7), C64::ZERO, C64::new(-5e-7,1e-7)]
}
fn near(a: C64, b: C64) {
    assert!((a-b).abs() <= 3e-8*b.abs().max(1e-23), "{a:?} != {b:?}");
}
fn near_power(a: f64, b: f64) {
    assert!((a-b).abs() <= 3e-8*b.abs().max(1e-26), "{a} != {b}");
}

#[test]
fn full_loaded_shell_matches_original_force_and_independent_combined_boundary_source() {
    with_cx(|cx,_| {
        let shell = source(mesh(),cx);
        assert_eq!(shell.model().free,6);
        for i in 0..6 { assert_eq!(shell.model().dof_map[18+i],Some(i)); }
        let geometry = TensorShellRadiation::new(&shell,medium(),geometry_controls(),cx).unwrap();
        assert_eq!(geometry.surface().areas().len(),16);
        assert_eq!(geometry.normal_velocity_rows().len(),6);
        assert!(geometry.normal_velocity_rows()[3..].iter().flatten().all(|x| *x==0.0));
        // Even the completely clamped bottom facet remains in the acoustic surface.
        assert_eq!(geometry.surface().centroids().iter().filter(|p| p[2]==0.0).count(),4);
        let prepared = geometry.prepare(200.0,boundary_budget(),cx).unwrap();
        let loaded = prepared.prepare_frequency(100.0,frequency_controls(),cx).unwrap();
        let result = loaded.solve(&force(),cx).unwrap();
        let motion = result.motion();
        assert_eq!(motion.displacement.len(),6);
        assert!(motion.average_viscous_power_w > 0.0);
        assert!(result.power().outgoing_power_w > 0.0);
        assert!(result.relative_coupled_power_defect() < 1e-8);
        assert!(result.relative_port_work_defect() < 1e-8);
        assert!(result.power().relative_flux_defect < 1e-8);
        near_power(motion.average_input_power_w,
            motion.average_viscous_power_w + result.power().outgoing_power_w);

        // Directly call the native full material/physical-port factor owner.
        let mut identity = vec![0.0;36];
        for i in 0..6 { identity[6*i+i]=1.0; }
        let direct_shell = shell.prepare_harmonic_response_with_port_impedance(
            100.0,&identity,6,loaded.impedance(),frequency_controls().mechanical,
            frequency_controls().ports,cx).unwrap().solve(&force(),cx).unwrap();
        for (&a,&b) in motion.displacement.iter().zip(&direct_shell.displacement) { near(a,b); }
        let bare = shell.prepare_harmonic_response(100.0,frequency_controls().mechanical,cx)
            .unwrap().solve(&force(),cx).unwrap();
        let bare_scale = bare.displacement.iter().map(|x| x.abs()).fold(0.0,f64::max);
        assert!(motion.displacement.iter().zip(&bare.displacement)
            .any(|(&a,&b)| (a-b).abs() > 1e-4*bare_scale));

        // Independent P1 interpolation: only vertex3 moves, and its shape on
        // every tetrahedron face containing it is z/0.5. Do not read adapter rows.
        let velocity: Vec<_> = geometry.surface().centroids().iter()
            .zip(geometry.surface().normals()).map(|(point,normal)| {
                (0..3).fold(C64::ZERO,|sum,d| sum+motion.velocity[d].scale(normal[d]))
                    .scale(point[2]/0.5)
            }).collect();
        for (&a,&b) in result.field().surface_velocity_m_s().iter().zip(&velocity) { near(a,b); }
        let direct = prepared.boundary().solve_batch(loaded.wavenumber_per_m(),medium(),
            &[&velocity],||false).unwrap();
        for (&a,&b) in result.field().surface_pressure_pa().iter().zip(&direct[0].pressure) { near(a,b); }
        near(result.power().surface_power_w,direct[0].surface_power);
        near_power(result.power().outgoing_power_w,direct[0].outgoing_power_w);
        let observers = loaded.prepare_receivers(&[[0.1,-0.3,0.8],[-0.3,0.1,0.2]],
            0.01,Options::default(),true,cx).unwrap();
        let observed = observers.evaluate_coherent_velocity(result.field()).unwrap();
        let reference = observers.evaluate_velocity(&direct[0]).unwrap();
        for r in 0..2 {
            near(observed.scalar.pressure[r],reference.scalar.pressure[r]);
            for d in 0..3 { near(observed.particle_velocity_m_s[r][d],reference.particle_velocity_m_s[r][d]); }
        }
    });
}

#[test]
fn repeated_forces_share_loaded_sources_and_coherent_directional_arrays() {
    with_cx(|cx,_| {
        let shell = source(mesh(),cx);
        let geometry = TensorShellRadiation::new(&shell,medium(),geometry_controls(),cx).unwrap();
        let prepared = geometry.prepare(200.0,boundary_budget(),cx).unwrap();
        // Receiver geometry is admitted before preparing any acoustic source column.
        let receiver_geometry = prepared.receiver_geometry(&[[0.1,-0.3,0.8],[0.1,-0.3,0.8]],
            0.01,cx).unwrap();
        let loaded = prepared.prepare_frequency(90.0,frequency_controls(),cx).unwrap();
        let z = loaded.impedance().to_vec();
        let observer = receiver_geometry.prepare(loaded.wavenumber_per_m(),medium(),
            Options::default(),true,||false).unwrap();
        let gain = C64::new(-0.7,0.4);
        let first = loaded.solve(&force(),cx).unwrap();
        let scaled_force = force().map(|x| x*gain);
        let second = loaded.solve(&scaled_force,cx).unwrap();
        for (&a,&b) in second.motion().displacement.iter().zip(&first.motion().displacement) {
            near(a,b*gain);
        }
        near_power(second.power().outgoing_power_w,first.power().outgoing_power_w*gain.norm_sq());
        near_power(second.motion().average_viscous_power_w,
            first.motion().average_viscous_power_w*gain.norm_sq());
        let work = observer.kernel_evaluations();
        let a = observer.evaluate_coherent_velocity(first.field()).unwrap();
        let b = observer.evaluate_coherent_velocity(second.field()).unwrap();
        let front = FirstOrder::new(0.5,[0.0,0.0,1.0]).unwrap();
        let back = FirstOrder::new(0.5,[0.0,0.0,-1.0]).unwrap();
        for r in 0..2 {
            near(b.scalar.pressure[r],a.scalar.pressure[r]*gain);
            for d in 0..3 { near(b.particle_velocity_m_s[r][d],a.particle_velocity_m_s[r][d]*gain); }
            near(front.observe(a.scalar.pressure[r],a.particle_velocity_m_s[r],medium()).unwrap()
                +back.observe(a.scalar.pressure[r],a.particle_velocity_m_s[r],medium()).unwrap(),
                a.scalar.pressure[r]);
        }
        assert_eq!(work,observer.kernel_evaluations());
        assert_eq!(loaded.impedance(),z.as_slice());
        let zero = loaded.solve(&[C64::ZERO;6],cx).unwrap();
        assert_eq!(zero.power().outgoing_power_w,0.0);
        assert_eq!(zero.power().surface_power_w,C64::ZERO);
        assert_eq!(zero.relative_coupled_power_defect(),0.0);
        assert_eq!(zero.relative_port_work_defect(),0.0);
    });
}

#[test]
fn original_geometry_material_bands_complete_coordinates_and_resources_refuse() {
    with_cx(|cx,gate| {
        let shell = source(mesh(),cx);
        let original_mesh = shell.mesh().clone();
        assert!(matches!(TensorShellRadiation::new(&shell,medium(),
            TensorShellRadiationControls { maximum_coordinates:5, ..geometry_controls() },cx),
            Err(TensorShellRadiationError::Budget { .. })));
        assert!(TensorShellRadiation::new(&shell,medium(),TensorShellRadiationControls {
            maximum_setup_bytes:1, ..geometry_controls()
        },cx).is_err());
        assert!(TensorShellRadiation::new(&shell,medium(),TensorShellRadiationControls {
            maximum_panels:15, ..geometry_controls()
        },cx).is_err());
        let geometry = TensorShellRadiation::new(&shell,medium(),geometry_controls(),cx).unwrap();
        let prepared = geometry.prepare(200.0,boundary_budget(),cx).unwrap();
        assert!(prepared.receiver_geometry(&[],0.01,cx).is_err());
        assert!(prepared.receiver_geometry(&[[0.1,0.1,0.1]],0.01,cx).is_err());
        assert!(prepared.prepare_frequency(201.0,frequency_controls(),cx).is_err());
        assert!(prepared.prepare_frequency(100.0,TensorShellRadiationFrequencyControls {
            maximum_boundary_payload_bytes:1, ..frequency_controls()
        },cx).is_err());
        assert!(prepared.prepare_frequency(100.0,TensorShellRadiationFrequencyControls {
            field:CoherentPortFieldLimits { maximum_payload_bytes:1, ..frequency_controls().field },
            ..frequency_controls()
        },cx).is_err());
        assert!(matches!(prepared.prepare_frequency(180.0,frequency_controls(),cx),
            Err(TensorShellRadiationError::Shell(TensorShellHarmonicError::Material { .. }))));
        let loaded = prepared.prepare_frequency(100.0,frequency_controls(),cx).unwrap();
        assert!(matches!(loaded.solve(&force().map(|x| x.scale(1e12)),cx),
            Err(TensorShellRadiationError::Shell(TensorShellHarmonicError::Material { .. }))));
        assert_eq!(shell.mesh(),&original_mesh);

        let mut open = mesh();
        open.tris.pop();
        let open_source = source(open,cx);
        let open_geometry = TensorShellRadiation::new(&open_source,medium(),geometry_controls(),cx).unwrap();
        assert!(open_geometry.prepare(200.0,boundary_budget(),cx).is_err());
        let mut inward = mesh();
        for triangle in &mut inward.tris { triangle.swap(1,2); }
        let inward_source = source(inward,cx);
        let inward_geometry = TensorShellRadiation::new(&inward_source,medium(),geometry_controls(),cx).unwrap();
        assert!(inward_geometry.prepare(200.0,boundary_budget(),cx).is_err());
        gate.request();
        assert!(matches!(loaded.solve(&force(),cx),Err(TensorShellRadiationError::Cancelled)));
    });
}
