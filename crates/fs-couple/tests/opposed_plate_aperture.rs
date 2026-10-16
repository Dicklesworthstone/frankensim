//! Opposed-blade mechanics and shared-owner integration, not measured-oboe claims.
use fs_couple::bernoulli_aperture::dynamic::{ApertureDrive, ApertureState, DynamicAperture};
use fs_couple::bernoulli_aperture::dynamic::force::PlateForceFootprint;
use fs_couple::bernoulli_aperture::dynamic::relaxation::{
    InitialApertureMemory, PlateRelaxationRegion, PlateRelaxationSpec,
};
use fs_couple::bernoulli_aperture::plate::{
    PlateApertureOptions, PlateApertureReduction, PlateApertureSymmetry,
};
use fs_couple::bernoulli_aperture::plate::closure::PlateClosureSpec;
use fs_exec::CancelGate;
use fs_material::visco::GeneralizedMaxwell;
use fs_plate::{AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateSection, SliceOptions};

fn blade(e: f64, density: f64) -> PlateApertureReduction {
    let mesh = PlateMesh::rectangle(0.025, 0.01, 6, 2);
    let root = mesh.nodes.iter().enumerate().filter(|(_, p)| p.0 == 0.0)
        .map(|(i, _)| i).collect();
    let edges = mesh.boundary_edges().into_iter().filter(|&(a, b)|
        (mesh.nodes[a].0 - 0.025).abs() < 1e-12 && (mesh.nodes[b].0 - 0.025).abs() < 1e-12)
        .map(|(a, b)| [a, b]).collect();
    let section = PlateSection::isotropic(e, 0.3, 0.0003, density).unwrap();
    let chart = PlateChart::with_boundary_and_regions(mesh, section, root, vec![]).unwrap();
    PlateApertureReduction::from_chart(chart, PlateApertureOptions {
        assembly: AssemblyOptions { pretension: 0.0, support: EdgeSupport::Clamped },
        eigenvalue_window: (1.0, (core::f64::consts::TAU * 450.0).powi(2)), mode_index: 0,
        eigensolver: SliceOptions::default(), max_nodes: 100, max_triangles: 100,
        slit_edges: edges, rest_opening_m: 0.0002, damping_ratio: 0.0,
        max_slit_mode_variation: 0.25, max_slope: 0.1,
    }, &CancelGate::new()).unwrap()
}
fn near(actual: f64, expected: f64, tol: f64) {
    assert!((actual - expected).abs() <= tol * expected.abs().max(1e-30),
        "{actual:e} != {expected:e}");
}
fn closure(r: &PlateApertureReduction) -> PlateClosureSpec {
    PlateClosureSpec {
        nodal_rest_gap_m: vec![r.options().rest_opening_m; r.chart().mesh.node_count()],
        lay_triangles: (0..r.chart().mesh.tris.len()).collect(),
        stiffness_pa_per_m_alpha: 1e10, alpha: 2.0, internal_loss_s_per_m: 0.0,
        provenance: "authored common-interface comparison, not a measured material pair".into(),
        max_penetration_m: 0.001,
    }
}
fn with_memory(r: PlateApertureReduction) -> DynamicAperture {
    let count = r.chart().mesh.tris.len();
    let state = ApertureState {
        opening_m: r.options().rest_opening_m - 1e-5, opening_velocity_m_s: 0.0,
    };
    let contact = closure(&r);
    DynamicAperture::from_plate_with_closure(r, contact, 1.2, 1e6, 1e-5, 128, state).unwrap()
        .with_plate_relaxation(PlateRelaxationSpec {
            regions: vec![PlateRelaxationRegion {
                triangles: (0..count).collect(),
                material: GeneralizedMaxwell::new(4e9, vec![(2e9, 0.001)]).unwrap(),
                poisson_ratio: 0.3, band_hz: (0.0, 10000.0),
                provenance: "authored proportional material for symmetry regression".into(),
            }], max_branches: 1, max_dt_over_tau: 0.1, max_angular_step: 0.5,
        }, InitialApertureMemory::Unrelaxed).unwrap()
}

#[test]
fn opposed_reduction_is_the_sum_of_two_physical_plate_quadratic_forms() {
    let single = blade(4e9, 900.0);
    let pair = single.clone().into_symmetric_opposed_pair().unwrap();
    assert_eq!(pair.symmetry(), PlateApertureSymmetry::SymmetricOpposedPair);
    assert_eq!(single.moving_blade_count(), 1);
    assert_eq!(pair.moving_blade_count(), 2);
    assert!(pair.clone().into_symmetric_opposed_pair().is_err());
    let model = pair.chart().assemble(&[], &pair.options().assembly).unwrap();
    let mut shape = vec![0.0; model.free];
    for (node, value) in pair.shape_per_opening().iter().enumerate() {
        for c in 0..3 { if let Some(i) = model.dof_map[3*node+c] { shape[i] = value[c]; } }
    }
    let mut product = vec![0.0; shape.len()];
    model.m.spmv(&shape, &mut product);
    near(pair.mass_kg(), 2.0 * shape.iter().zip(&product).map(|(a,b)| a*b).sum::<f64>(), 1e-10);
    model.k.spmv(&shape, &mut product);
    near(pair.stiffness_n_m(), 2.0 * shape.iter().zip(&product).map(|(a,b)| a*b).sum::<f64>(), 1e-6);
    let mut per_blade_area = 0.0;
    for tri in &pair.chart().mesh.tris {
        let [a,b,c] = tri.map(|i| pair.chart().mesh.nodes[i]);
        let da = 0.5 * ((b.0-a.0)*(c.1-a.1)-(b.1-a.1)*(c.0-a.0));
        per_blade_area += da * tri.iter().map(|&i| pair.shape_per_opening()[i][0]).sum::<f64>() / 3.0;
    }
    near(pair.pressure_area_m2(), 2.0 * per_blade_area, 1e-12);
    assert_eq!(pair.mass_kg(), 0.5 * single.mass_kg());
    assert_eq!(pair.stiffness_n_m(), 0.5 * single.stiffness_n_m());
    assert_eq!(pair.closing_pressure_pa(), 0.5 * single.closing_pressure_pa());
    assert_eq!(pair.pressure_area_m2(), single.pressure_area_m2());
    assert_eq!(pair.width_m(), single.width_m());
    assert_eq!(pair.options().rest_opening_m, single.options().rest_opening_m);
    assert_eq!(pair.eigenvalue_interval(), single.eigenvalue_interval());
    near(pair.stiffness_n_m()/pair.mass_kg(), single.stiffness_n_m()/single.mass_kg(), 1e-14);
}

#[test]
fn opposed_closure_uses_total_gap_and_counts_the_common_interface_once() {
    let single = blade(4e9, 900.0);
    let pair = single.clone().into_symmetric_opposed_pair().unwrap();
    let mut profile = closure(&single);
    // Nonuniform slit crosses zero partway along an edge as the blades close.
    for (i, gap) in profile.nodal_rest_gap_m.iter_mut().enumerate() {
        *gap *= 0.7 + 0.6 * single.chart().mesh.nodes[i].1 / 0.01;
    }
    let a = single.compile_closure(profile.clone()).unwrap();
    let b = pair.compile_closure(profile).unwrap();
    assert_eq!(a.contact_law().weights(), b.contact_law().weights());
    assert_eq!(a.contact_nodes(), b.contact_nodes());
    let h = single.options().rest_opening_m;
    for opening in [h, 0.5*h, 0.0, -0.2*h] {
        assert_eq!(a.probe(opening).unwrap(), b.probe(opening).unwrap());
        for node in 0..single.chart().mesh.node_count() {
            assert_eq!(a.nodal_gap_m(node, opening).unwrap(), b.nodal_gap_m(node, opening).unwrap());
        }
    }
}

#[test]
fn physical_leaf_motion_and_equally_split_total_force_are_work_conjugate() {
    let single = blade(4e9, 900.0);
    let pair = single.clone().into_symmetric_opposed_pair().unwrap();
    let node = single.options().slit_edges[0][0];
    let opening = single.options().rest_opening_m - 2e-5;
    let speed = 0.02;
    let a = single.nodal_motion(node, opening, speed).unwrap();
    let b = pair.nodal_motion(node, opening, speed).unwrap();
    for i in 0..3 {
        assert_eq!(b.displacement_rotation[i], 0.5*a.displacement_rotation[i]);
        assert_eq!(b.velocity_rotation_rate[i], 0.5*a.velocity_rotation_rate[i]);
    }
    assert_eq!(pair.max_slope_at(opening), 0.5*single.max_slope_at(opening));
    let port = pair.force_port(PlateForceFootprint::Node(node)).unwrap();
    let total_force = 0.8;
    near(port.generalized_force_n(total_force).unwrap()*speed,
        2.0*(total_force/2.0)*b.velocity_rotation_rate[0], 1e-14);
    near(port.velocity_m_s(speed).unwrap(), b.velocity_rotation_rate[0], 1e-14);
    let patch = pair.force_port(PlateForceFootprint::Patch((0..pair.chart().mesh.tris.len()).collect())).unwrap();
    // Two equal face pressures recover the SAME combined pressure-area port.
    let pressure = 42.0;
    near(patch.generalized_force_n(2.0*patch.area_m2().unwrap()*pressure).unwrap(),
        pressure*pair.pressure_area_m2(), 1e-12);
}

#[test]
fn opposed_material_memory_counts_both_blades_and_preserves_the_loaded_energy_balance() {
    let single = blade(4e9, 900.0);
    let pair = single.clone().into_symmetric_opposed_pair().unwrap();
    let a = with_memory(single);
    let mut b = with_memory(pair);
    let am = a.relaxation().unwrap();
    let bm = b.relaxation().unwrap();
    near(bm.region_bending_stiffness_n_m()[0], 0.5*am.region_bending_stiffness_n_m()[0], 1e-12);
    near(bm.branches()[0].stiffness, 0.5*am.branches()[0].stiffness, 1e-12);
    near(bm.stored_energy_j(), 0.5*am.stored_energy_j(), 1e-12);
    assert_eq!(bm.branches()[0].relaxation_time_s, am.branches()[0].relaxation_time_s);
    near(b.stored_energy_j(), 0.5*a.stored_energy_j(), 1e-12);
    let area = b.plate_reduction().unwrap().pressure_area_m2();
    for i in 0..64 {
        let old = b.state();
        let frame = b.step(ApertureDrive {
            upstream_pressure_pa: 10.0 + 0.01*i as f64,
            incoming_pressure_pa: 0.0, body_flow_m3_s: 0.0,
        }).unwrap();
        let vm = 0.5*(old.opening_velocity_m_s+frame.state.opening_velocity_m_s);
        near(frame.swept_flow_m3_s, -area*vm, 1e-12);
        assert!(frame.flow_residual_m3_s.abs() < 1e-11);
        assert!(frame.balance_residual_j().abs() < 1e-11);
        assert!(frame.relaxation_loss_j >= 0.0);
        assert_eq!(b.plate_reduction().unwrap().moving_blade_count(), 2);
    }
    let state = b.state();
    let memory = b.relaxation().unwrap().memory_sqrt_j().to_vec();
    assert!(b.step(ApertureDrive { upstream_pressure_pa: f64::NAN,
        incoming_pressure_pa: 0.0, body_flow_m3_s: 0.0 }).is_err());
    assert_eq!(b.state(), state);
    assert_eq!(b.relaxation().unwrap().memory_sqrt_j(), memory.as_slice());
    assert_eq!(b.accepted_steps(), 64);
}

#[test]
fn opposed_material_swaps_change_elasticity_and_inertia_without_retuning_geometry() {
    let a = blade(4e9, 900.0).into_symmetric_opposed_pair().unwrap();
    let stiff = blade(8e9, 900.0).into_symmetric_opposed_pair().unwrap();
    let heavy = blade(4e9, 1800.0).into_symmetric_opposed_pair().unwrap();
    near(stiff.stiffness_n_m()/a.stiffness_n_m(), 2.0, 2e-5);
    near(stiff.mass_kg()/a.mass_kg(), 1.0, 2e-5);
    near(heavy.mass_kg()/a.mass_kg(), 2.0, 2e-5);
    near(heavy.stiffness_n_m()/a.stiffness_n_m(), 1.0, 2e-5);
    near(stiff.closing_pressure_pa()/a.closing_pressure_pa(), 2.0, 2e-5);
    near(heavy.closing_pressure_pa()/a.closing_pressure_pa(), 1.0, 2e-5);
    assert_eq!(stiff.width_m(), a.width_m());
    assert_eq!(heavy.options().rest_opening_m, a.options().rest_opening_m);
}
