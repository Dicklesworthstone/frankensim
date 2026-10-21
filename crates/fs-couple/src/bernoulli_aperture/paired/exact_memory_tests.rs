//! Synthetic exact-material coupling checks, not an oboe fidelity certificate.
use super::*;
use crate::bernoulli_aperture::dynamic::relaxation::PlateRelaxationRegion;
use crate::bernoulli_aperture::plate::PlateApertureOptions;
use fs_exec::CancelGate;
use fs_material::{gas::GasSpec, visco::GeneralizedMaxwell};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateSection, SliceOptions};

const DT: f64 = 1e-5;
fn plate(e: f64, density: f64) -> PlateApertureReduction {
    let mesh = PlateMesh::rectangle(0.025, 0.01, 6, 2);
    let root = mesh.nodes.iter().enumerate().filter(|(_, p)| p.0 == 0.0).map(|(i, _)| i).collect();
    let edges = mesh.boundary_edges().into_iter().filter(|&(a, b)|
        (mesh.nodes[a].0-0.025).abs() < 1e-12 && (mesh.nodes[b].0-0.025).abs() < 1e-12)
        .map(|(a, b)| [a, b]).collect();
    let section = PlateSection::isotropic(e, 0.3, 0.0003, density).unwrap();
    let chart = PlateChart::with_boundary_and_regions(mesh, section, root, vec![]).unwrap();
    PlateApertureReduction::from_chart(chart, PlateApertureOptions {
        assembly: AssemblyOptions { pretension: 0.0, support: EdgeSupport::Clamped },
        eigenvalue_window: (1.0, (core::f64::consts::TAU*450.0).powi(2)), mode_index: 0,
        eigensolver: SliceOptions::default(), max_nodes: 100, max_triangles: 100,
        slit_edges: edges, rest_opening_m: 0.0002, damping_ratio: 0.0,
        max_slit_mode_variation: 0.25, max_slope: 0.1,
    }, &CancelGate::new_clock_free()).unwrap()
}
fn build(terms: Option<[Vec<(f64, f64)>; 2]>) -> PairedAperture {
    let plates = [plate(4e9, 900.0), plate(8e9, 1300.0)];
    let laws = terms.map(|terms| std::array::from_fn::<_, 2, _>(|i| PlateRelaxationSpec {
        regions: vec![PlateRelaxationRegion {
            triangles: (0..plates[i].chart().mesh.tris.len()).collect(),
            material: GeneralizedMaxwell::new([4e9, 8e9][i], terms[i].clone()).unwrap(),
            poisson_ratio: 0.3, band_hz: (0.0, 1000.0),
            provenance: "synthetic independent exact blade spectrum".into(),
        }], max_branches: 64, max_dt_over_tau: 0.1, max_angular_step: 0.1,
    }));
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let contact = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 0.0, 2.0,
        "explicit no-contact exact-material comparison".into()).unwrap();
    let mut model = PairedAperture::from_plates_with_propagator(plates, gas, 1e6, DT, 4096,
        PairedApertureLimits { maximum_abs_displacement_m: 0.001, maximum_abs_velocity_m_s: 10.0,
            maximum_energy_j: 1.0, max_evaluations: 256, velocity_tolerance_m_s: 1e-15,
            relative_tolerance: 1e-13 },
        [BladeState { displacement_m: -3e-6, velocity_m_s: 0.0 },
         BladeState { displacement_m: 2e-6, velocity_m_s: 0.0 }], contact,
        BladePropagator::ExactZoh).unwrap();
    if let Some([first, second]) = laws {
        model = model.with_blade_relaxation(0, first, InitialApertureMemory::Unrelaxed).unwrap()
            .with_blade_relaxation(1, second, InitialApertureMemory::Relaxed).unwrap();
    }
    model
}
fn material_pair() -> PairedAperture {
    build(Some([vec![(2e9, 0.001), (1e9, 0.01)], vec![(3e9, 0.003)]]))
}
fn near(a: f64, b: f64, scale: f64) {
    assert!((a-b).abs() <= 3e-8*scale.max(1e-20), "{a:e} != {b:e}");
}

#[test]
fn g1_independent_exact_memories_match_a_physical_momentum_state_generator() {
    let mut model = material_pair();
    let spec = *model.spec();
    // Independent coordinates [x, momentum, z..., integral(material force), F].
    // This does not reuse ExactBlade's scaled generator, endpoint map or loss forms.
    let maps: Vec<_> = (0..2).map(|i| {
        let b = spec.blades[i];
        let arms = model.blade_relaxation(i).unwrap().branches();
        let n = 4+arms.len();
        let mut a = vec![0.0; n*n];
        a[1] = DT/b.mass_kg;
        a[n] = -DT*b.stiffness_n_m;
        a[n+n-1] = DT;
        for (j, arm) in arms.iter().enumerate() {
            let c = arm.stiffness.sqrt()*arm.projection[0];
            a[n] -= DT*c*c;
            a[n+2+j] = DT*c;
            a[(2+j)*n] = DT*c/arm.relaxation_time_s;
            a[(2+j)*n+2+j] = -DT/arm.relaxation_time_s;
            a[(n-2)*n] += DT*c*c;
            a[(n-2)*n+2+j] = -DT*c;
        }
        fs_phs::linear_exact::exponential(&a, n).unwrap()
    }).collect();
    let mut states: Vec<Vec<f64>> = (0..2).map(|i| {
        let mut x = vec![model.blades[i].displacement_m,
            spec.blades[i].mass_kg*model.blades[i].velocity_m_s];
        x.extend_from_slice(model.blade_relaxation(i).unwrap().memory_sqrt_j());
        x.extend_from_slice(&[0.0, 0.0]);
        x
    }).collect();
    let energy = |i: usize, x: &[f64], arms: &[fs_phs::RelaxationBranch]| {
        let b = spec.blades[i];
        0.5*b.stiffness_n_m*x[0]*x[0]+0.5*x[1]*x[1]/b.mass_kg
            + arms.iter().enumerate().map(|(j, arm)|
                0.5*(arm.stiffness.sqrt()*arm.projection[0]*x[0]-x[2+j]).powi(2)).sum::<f64>()
    };
    for tick in 0..64 {
        let dp = if tick < 32 { 2.0_f64 } else { -1.0 };
        let forces = if tick < 40 { [1e-4, -2e-4] } else { [-1e-4, 3e-4] };
        let expected: Vec<Vec<f64>> = (0..2).map(|i| {
            let n = states[i].len();
            states[i][n-2] = 0.0;
            states[i][n-1] = forces[i]-spec.blades[i].pressure_area_m2*dp;
            (0..n).map(|r| (0..n).map(|c| maps[i][r*n+c]*states[i][c]).sum()).collect()
        }).collect();
        let mean_v: [f64; 2] = std::array::from_fn(|i| (expected[i][0]-states[i][0])/DT);
        let gap_before = spec.rest_opening_m+states[0][0]+states[1][0];
        let gap_after = spec.rest_opening_m+expected[0][0]+expected[1][0];
        let jet = spec.width_m*f64::midpoint(gap_before, gap_after).max(0.0)*dp.signum()
            *(2.0*dp.abs()/spec.density_kg_m3).sqrt();
        let swept = -(0..2).map(|i| spec.blades[i].pressure_area_m2*mean_v[i]).sum::<f64>();
        let mut loss = 0.0;
        for i in 0..2 {
            let m = model.blade_relaxation(i).unwrap();
            let n = states[i].len();
            let held = states[i][n-1];
            near(model.exact.as_ref().unwrap()[i].mean_memory_force_n(
                model.blades[i], held, DT, m.memory_sqrt_j()), expected[i][n-2]/DT, 0.01);
            loss += held*(expected[i][0]-states[i][0])
                - energy(i, &expected[i], m.branches())+energy(i, &states[i], m.branches());
        }
        let frame = model.step_with_forces(ApertureDrive {
            upstream_pressure_pa: dp+0.6+spec.impedance_pa_s_m3*(jet+swept+1e-8),
            incoming_pressure_pa: 0.3, body_flow_m3_s: 1e-8,
        }, forces).unwrap();
        near(frame.pressure_drop_pa, dp, dp.abs());
        for i in 0..2 {
            near(frame.blades[i].displacement_m, expected[i][0], spec.rest_opening_m);
            near(frame.blades[i].velocity_m_s, expected[i][1]/spec.blades[i].mass_kg, 0.1);
            assert!(frame.momentum_residual_n[i].abs() < 1e-9);
            let m = model.blade_relaxation(i).unwrap();
            for (j, z) in m.memory_sqrt_j().iter().enumerate() { near(*z, expected[i][2+j], 1e-3); }
        }
        near(frame.aperture.relaxation_loss_j, loss, frame.aperture.stored_energy_j+loss.abs());
        assert!(frame.aperture.relaxation_loss_j >= 0.0);
        let scale = frame.aperture.stored_energy_j+frame.aperture.storage_change_j.abs()
            +frame.aperture.dissipated_energy_j+frame.aperture.pressure_work_j.abs()
            +frame.aperture.mechanical_work_j.abs();
        near(frame.aperture.balance_residual_j(), 0.0, scale);
        states = expected;
    }
}

#[test]
fn g0_empty_exact_spectra_preserve_the_elastic_trajectory() {
    let mut a = build(None);
    let mut b = build(Some([vec![], vec![]]));
    for tick in 0..128 {
        let d = ApertureDrive { upstream_pressure_pa: if tick < 64 { 3.0 } else { -1.0 },
            ..ApertureDrive::default() };
        assert_eq!(a.step(d).unwrap(), b.step(d).unwrap());
    }
}

#[test]
fn g3_changing_only_one_blades_relaxation_time_changes_coupled_pressure() {
    let mut a = build(Some([vec![], vec![(3e9, 0.001)]]));
    let mut b = build(Some([vec![], vec![(3e9, 0.02)]]));
    assert_eq!(a.spec(), b.spec());
    let mut changed = false;
    for tick in 0..256 {
        let drive = ApertureDrive { upstream_pressure_pa: if tick < 128 { 5.0 } else { 0.0 },
            ..ApertureDrive::default() };
        let x = a.step(drive).unwrap();
        let y = b.step(drive).unwrap();
        changed |= (x.aperture.bore_pressure_pa-y.aperture.bore_pressure_pa).abs() > 1e-10;
    }
    assert!(changed);
    assert!(a.blade_relaxation(1).unwrap().dissipated_energy_j() > 0.0);
    assert!(b.blade_relaxation(1).unwrap().dissipated_energy_j() > 0.0);
}

#[test]
fn g0_exact_trial_and_refusal_do_not_publish_either_blades_history() {
    let mut model = material_pair();
    let before = *model.blade_states();
    let energy = model.stored_energy_j();
    let histories: Vec<_> = (0..2).map(|i|
        model.blade_relaxation(i).unwrap().memory_sqrt_j().to_vec()).collect();
    let drive = ApertureDrive { upstream_pressure_pa: 3.0, ..ApertureDrive::default() };
    let trial = model.preview(drive, [1e-4, -1e-4]).unwrap();
    assert!(model.step_with_forces(drive, [f64::NAN, 0.0]).is_err());
    assert_eq!(model.blade_states(), &before);
    assert_eq!(model.stored_energy_j(), energy);
    assert_eq!(model.accepted_steps(), 0);
    for i in 0..2 {
        assert_eq!(model.blade_relaxation(i).unwrap().memory_sqrt_j(), histories[i]);
        assert_eq!(model.blade_relaxation(i).unwrap().dissipated_energy_j(), 0.0);
    }
    let expected = trial.frame;
    model.accept(trial);
    assert_eq!(model.stored_energy_j(), expected.aperture.stored_energy_j);
    assert_eq!(model.accepted_steps(), 1);
    let histories: Vec<_> = (0..2).map(|i|
        model.blade_relaxation(i).unwrap().memory_sqrt_j().to_vec()).collect();
    model.spec.max_steps = 1;
    assert!(model.step(drive).is_err());
    for i in 0..2 { assert_eq!(model.blade_relaxation(i).unwrap().memory_sqrt_j(), histories[i]); }
    model.extend_step_budget(2).unwrap();
    model.step(drive).unwrap();
}
