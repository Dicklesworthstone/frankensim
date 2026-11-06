use super::*;
use crate::{Duct, Segment};
use fs_material::gas::GasSpec;

const DT: f64 = 1.0 / 96_000.0;
const TOL: f64 = 1e-11;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

/// A viscothermal cylinder bore driven at its inlet port, closed at the far end. DOFs are
/// `2 E (p + 1) + 2 E Q P`.
fn bore(radius_m: f64, elements: usize, order: usize, poles: usize, quadrature: usize, cap: usize) -> DgDuct {
    let geometry = Duct {
        segments: vec![Segment::Cylinder { radius: radius_m, length: 0.57 }],
    };
    DgDuct::new(&geometry, &gas(), DgMesh { elements, order }, DgBoundary::Port(PortId(0)), DgBoundary::Closed)
        .unwrap()
        .with_dof_cap(cap)
        .unwrap()
        .with_spatial_viscothermal(DgWallResolution {
            poles_per_field: poles,
            quadrature_points: quadrature,
            band_rad_s: [core::f64::consts::TAU * 50.0, core::f64::consts::TAU * 4800.0],
            maximum_relative_error: 0.01,
        })
        .unwrap()
}

fn norm(v: &[f64]) -> f64 {
    det::sqrt(v.iter().map(|x| x * x).sum())
}

/// A held inlet flow: a 24-step raised-cosine pulse of 1e-6 m^3/s, then zero (never zero while it lasts).
fn effort(step: usize) -> f64 {
    if step < 24 {
        let x = core::f64::consts::PI * (step as f64 + 1.0) / 25.0;
        1.0e-6 * x.sin() * x.sin()
    } else {
        0.0
    }
}

fn input(u: f64) -> [PortInput; 1] {
    [PortInput { id: PortId(0), effort: u }]
}

fn kernel_of(propagator: &DgPropagator) -> &KrylovKernel {
    match &propagator.kernel {
        PropagatorKernel::Krylov(k) => k,
        _ => panic!("expected the Krylov kernel"),
    }
}

/// The dense and the Krylov propagators of the same duct agree step by step: state, port
/// outputs, the multiport response and the energy ledger.
fn assert_matches_dense(duct: &DgDuct, steps: usize) {
    let dense = duct.propagator(DT, TOL).unwrap();
    let sparse = duct.build_krylov_propagator(DT, TOL).unwrap();
    assert!(dense.krylov_report().is_none());
    let mut sd = duct.zero_state();
    let mut sk = duct.zero_state();
    let mut worst_state = 0.0_f64;
    let mut worst_rate = 0.0_f64;
    let mut worst_resp = 0.0_f64;
    for k in 0..steps {
        let rd = dense.multiport_response(&sd).unwrap();
        let rk = sparse.multiport_response(&sk).unwrap();
        let scale = rd.history_pressure_pa[0].abs() + rd.impedance_pa_s_m3[0].abs() * 1.0e-6 + 1e-30;
        worst_resp = worst_resp
            .max((rd.history_pressure_pa[0] - rk.history_pressure_pa[0]).abs() / scale)
            .max((rd.impedance_pa_s_m3[0] - rk.impedance_pa_s_m3[0]).abs() / rd.impedance_pa_s_m3[0].abs().max(1e-30));
        let od = dense.step(&mut sd, &input(effort(k))).unwrap();
        let ok = sparse.step(&mut sk, &input(effort(k))).unwrap();
        let diff: Vec<f64> = sd.tilde.iter().zip(&sk.tilde).map(|(a, b)| a - b).collect();
        worst_state = worst_state.max(norm(&diff) / norm(&sd.tilde).max(1e-30));
        let (a, b) = (od.port_outputs[0].conjugate_rate, ok.port_outputs[0].conjugate_rate);
        worst_rate = worst_rate.max((a - b).abs() / a.abs().max(1e-12));
        let (ea, eb) = (od.ledger.stored_energy_j, ok.ledger.stored_energy_j);
        assert!((ea - eb).abs() <= 1e-8 * ea.abs().max(1e-300), "energy at step {k}: {ea} vs {eb}");
    }
    println!("KRYLOV-MEASURE dense vs krylov over {steps} steps: state {worst_state:e}, port mean pressure {worst_rate:e}, response {worst_resp:e}, krylov {:?}", sparse.krylov_report());
    assert!(worst_state < 1e-8, "state {worst_state}");
    assert!(worst_rate < 1e-7, "port mean pressure {worst_rate}");
    assert!(worst_resp < 1e-7, "response {worst_resp}");
}

#[test]
fn krylov_matches_the_dense_propagator_step_by_step() {
    // 3 elements, order 3, 4 wall points, 10 poles: 2*3*4 + 2*3*4*10 = 264 DOFs.
    let duct = bore(0.001, 3, 3, 10, 4, MAX_DG_PROPAGATOR_DOFS);
    assert_eq!(duct.total_dofs(), 264);
    assert_matches_dense(&duct, 300);
}

#[test]
fn krylov_matches_the_dense_propagator_on_a_stiff_narrow_bore() {
    // A 0.2 mm bore: the wall poles are far stiffer (decay rates ~1e7 per second), so the Arnoldi
    // dimension and the substep count are not those of the 1 mm bore. 5 elements, order 2.
    let duct = bore(0.0002, 5, 2, 5, 3, MAX_DG_PROPAGATOR_DOFS);
    assert_matches_dense(&duct, 240);
}

#[test]
fn above_the_dense_cap_the_opted_in_duct_takes_the_krylov_kernel() {
    // 21 elements, order 3, 4 wall points, 12 poles: 2*21*4 + 2*21*4*12 = 2184 DOFs.
    let duct = bore(0.001, 21, 3, 12, 4, MAX_DG_SPARSE_DOFS);
    assert_eq!(duct.total_dofs(), 2184);
    assert!(duct.total_dofs() > MAX_DG_PROPAGATOR_DOFS);
    let propagator = duct.propagator(DT, TOL).unwrap();
    let report = propagator.krylov_report().expect("above the dense cap the kernel is Krylov");
    assert_eq!(report.dofs, 2184);
    assert!(report.build_estimate <= TOL, "{report:?}");
    assert!(report.nonzeros < 2184 * 2184 / 10, "{report:?}");

    // The mean port pressure the response promises is the one the step delivers; the response and
    // the step share one solve; a clone shares it too.
    let mut state = duct.zero_state();
    let mut previous = f64::INFINITY;
    let mut peak = 0.0_f64;
    for k in 0..120 {
        let before = kernel_of(&propagator).report().solves;
        let response = propagator.port_response(&state, PortId(0)).unwrap();
        let multi = propagator.multiport_response(&state).unwrap();
        assert_eq!(multi.history_pressure_pa[0].to_bits(), response.history_pressure_pa.to_bits());
        assert_eq!(kernel_of(&propagator).report().solves, before + 1, "two reads of one state are one solve");
        let u = effort(k);
        let mut trial = state.clone();
        let out_trial = propagator.step(&mut trial, &input(u)).unwrap();
        assert_eq!(kernel_of(&propagator).report().solves, before + 1, "a trial step on a clone reuses the solve");
        let out = propagator.step(&mut state, &input(u)).unwrap();
        assert_eq!(out.port_outputs[0].conjugate_rate.to_bits(), out_trial.port_outputs[0].conjugate_rate.to_bits());
        let promised = response.history_pressure_pa + response.impedance_pa_s_m3 * u;
        let delivered = out.port_outputs[0].conjugate_rate;
        peak = peak.max(delivered.abs());
        assert!((promised - delivered).abs() <= 1e-9 * (promised.abs() + delivered.abs()).max(1e-30), "step {k}: {promised} vs {delivered}");
        assert!(response.impedance_pa_s_m3 >= 0.0);
        if k >= 24 {
            // unforced: a lossy duct only loses energy
            assert!(out.ledger.stored_energy_j <= previous * (1.0 + 1e-12) + 1e-300, "step {k}: {} after {previous}", out.ledger.stored_energy_j);
        }
        previous = out.ledger.stored_energy_j;
    }
    assert!(peak > 0.0);
    let report = propagator.krylov_report().unwrap();
    println!("KRYLOV-MEASURE above the dense cap: {report:?}");
    assert!(report.max_dimension < STATE_MAX_DIM, "{report:?}");
}

#[test]
fn the_dense_kernel_stays_for_ducts_at_or_under_the_dense_cap_even_when_opted_in() {
    let duct = bore(0.001, 3, 3, 10, 4, MAX_DG_SPARSE_DOFS);
    assert_eq!(duct.dof_cap(), MAX_DG_SPARSE_DOFS);
    assert!(duct.propagator(DT, TOL).unwrap().krylov_report().is_none());
}

#[test]
fn a_duct_that_did_not_opt_in_refuses_above_the_dense_cap_as_before() {
    let geometry = Duct {
        segments: vec![Segment::Cylinder { radius: 0.001, length: 0.57 }],
    };
    let plain = DgDuct::new(&geometry, &gas(), DgMesh { elements: 21, order: 3 }, DgBoundary::Port(PortId(0)), DgBoundary::Closed).unwrap();
    assert_eq!(plain.dof_cap(), MAX_DG_PROPAGATOR_DOFS);
    let refused = plain.with_spatial_viscothermal(DgWallResolution {
        poles_per_field: 12,
        quadrature_points: 4,
        band_rad_s: [core::f64::consts::TAU * 50.0, core::f64::consts::TAU * 4800.0],
        maximum_relative_error: 0.01,
    });
    assert!(refused.is_err(), "2184 DOFs without the opt-in must refuse");
    let duct = DgDuct::new(&geometry, &gas(), DgMesh { elements: 6, order: 3 }, DgBoundary::Port(PortId(0)), DgBoundary::Closed).unwrap();
    assert!(duct.clone().with_dof_cap(MAX_DG_PROPAGATOR_DOFS - 1).is_err());
    assert!(duct.clone().with_dof_cap(MAX_DG_SPARSE_DOFS + 1).is_err());
    assert!(duct.with_dof_cap(MAX_DG_SPARSE_DOFS).is_ok());
}

#[test]
fn an_unconverged_run_is_reported_and_a_solve_refuses() {
    let duct = bore(0.001, 3, 3, 10, 4, MAX_DG_PROPAGATOR_DOFS);
    let propagator = duct.build_krylov_propagator(DT, TOL).unwrap();
    let kernel = kernel_of(&propagator);
    let probe = probe_states(duct.total_dofs()).remove(0);
    // 500 steps at once with 5 vectors cannot be resolved
    let run = kernel.exp_action(&probe, 500.0 * DT, 5, 1e-12, false).unwrap();
    assert!(!run.converged && run.estimate > 1e-12, "{}", run.estimate);
    // the same operator at the build step converges well inside the budget
    let ok = kernel.exp_action(&probe, DT, STATE_MAX_DIM, STEP_TOL_FRACTION * TOL, false).unwrap();
    assert!(ok.converged);
}

#[test]
fn a_zero_state_solves_to_zero_and_a_nonfinite_state_refuses() {
    let duct = bore(0.001, 3, 3, 10, 4, MAX_DG_PROPAGATOR_DOFS);
    let propagator = duct.build_krylov_propagator(DT, TOL).unwrap();
    let kernel = kernel_of(&propagator);
    let n = duct.total_dofs();
    let (next, int, _) = kernel.solve_with(&vec![0.0; n], 1).unwrap();
    assert!(next.iter().chain(&int).all(|v| *v == 0.0));
    let mut bad = vec![0.0; n];
    bad[3] = f64::NAN;
    assert!(kernel.solve_with(&bad, 1).is_err());
}

#[test]
fn an_oboe_sized_bore_certifies_and_steps_on_the_krylov_kernel() {
    // The target configuration of the sparse kernel: a reed tube, a long cone and a bell, 39 elements, order 1,
    // 2 wall points, 40 poles: 2*39*2 + 2*39*2*40 = 6396 DOFs (three times the dense cap).
    let geometry = Duct {
        segments: vec![
            Segment::Cylinder { radius: 0.00087, length: 0.007 },
            Segment::Cone { inlet_radius: 0.00087, outlet_radius: 0.00825, length: 0.532 },
            Segment::Cone { inlet_radius: 0.00825, outlet_radius: 0.019, length: 0.095 },
        ],
    };
    let duct = DgDuct::new(&geometry, &gas(), DgMesh { elements: 39, order: 1 }, DgBoundary::Port(PortId(0)), DgBoundary::Closed)
        .unwrap()
        .with_dof_cap(MAX_DG_SPARSE_DOFS)
        .unwrap()
        .with_spatial_viscothermal(DgWallResolution {
            poles_per_field: 40,
            quadrature_points: 2,
            band_rad_s: [core::f64::consts::TAU * 50.0, core::f64::consts::TAU * 4800.0],
            maximum_relative_error: 0.01,
        })
        .unwrap();
    assert_eq!(duct.total_dofs(), 6396);
    let propagator = duct.propagator(DT, TOL).unwrap();
    let mut state = duct.zero_state();
    let mut previous = f64::INFINITY;
    for k in 0..40 {
        let response = propagator.port_response(&state, PortId(0)).unwrap();
        let u = effort(k);
        let out = propagator.step(&mut state, &input(u)).unwrap();
        let promised = response.history_pressure_pa + response.impedance_pa_s_m3 * u;
        let delivered = out.port_outputs[0].conjugate_rate;
        assert!((promised - delivered).abs() <= 1e-9 * (promised.abs() + delivered.abs()).max(1e-30), "step {k}: {promised} vs {delivered}");
        if k >= 24 {
            assert!(out.ledger.stored_energy_j <= previous * (1.0 + 1e-12) + 1e-300, "step {k}");
        }
        previous = out.ledger.stored_energy_j;
    }
    let report = propagator.krylov_report().expect("6396 DOFs take the Krylov kernel");
    println!("KRYLOV-MEASURE oboe-sized bore: {report:?}");
    assert_eq!(report.dofs, 6396);
    assert!(report.build_estimate <= TOL && report.max_dimension < STATE_MAX_DIM, "{report:?}");
}

fn wall_resolution(ceiling: usize, quadrature: usize) -> DgWallResolution {
    DgWallResolution {
        poles_per_field: ceiling,
        quadrature_points: quadrature,
        band_rad_s: [core::f64::consts::TAU * 50.0, core::f64::consts::TAU * 4800.0],
        maximum_relative_error: 0.01,
    }
}

/// A port-driven, closed-ended duct on `segments`, wall memory attached with a fixed pole count or with the
/// per-point minimum under that ceiling.
fn walled(segments: Vec<Segment>, elements: usize, quadrature: usize, ceiling: usize, min_poles: bool) -> DgDuct {
    let duct = DgDuct::new(&Duct { segments }, &gas(), DgMesh { elements, order: 1 }, DgBoundary::Port(PortId(0)), DgBoundary::Closed)
        .unwrap()
        .with_dof_cap(MAX_DG_SPARSE_DOFS)
        .unwrap();
    let resolution = wall_resolution(ceiling, quadrature);
    if min_poles {
        duct.with_spatial_viscothermal_min_poles(resolution).unwrap()
    } else {
        duct.with_spatial_viscothermal(resolution).unwrap()
    }
}

fn cone() -> Vec<Segment> {
    vec![Segment::Cone { inlet_radius: 0.00087, outlet_radius: 0.019, length: 0.634 }]
}

fn oboe_like() -> Vec<Segment> {
    vec![
        Segment::Cylinder { radius: 0.00087, length: 0.007 },
        Segment::Cone { inlet_radius: 0.00087, outlet_radius: 0.00825, length: 0.532 },
        Segment::Cone { inlet_radius: 0.00825, outlet_radius: 0.019, length: 0.095 },
    ]
}

/// The largest diagonal decay rate of the kernel's generator times dt: the stiffness the Arnoldi dimension follows.
fn stiffness(propagator: &DgPropagator) -> f64 {
    let g = &kernel_of(propagator).generator;
    let mut worst = 0.0_f64;
    for i in 0..g.n {
        for k in g.row_start[i]..g.row_start[i + 1] {
            if g.col[k] as usize == i {
                worst = worst.max(g.val[k].abs());
            }
        }
    }
    worst * DT
}

#[test]
fn min_poles_is_the_smallest_resolving_count_on_a_uniform_bore() {
    // Every point of a uniform bore has the same radius and gas, so the per-point minimum is one count: the first
    // count an independent scan from 1 finds. 3 elements, order 1, 4 wall points.
    let radius = 0.00087;
    let g = gas();
    let law = |poles: usize| {
        fs_phs::ZkDiffusive::new(radius, g.density_kg_m3(), g.dynamic_viscosity, g.gamma, g.prandtl, poles, wall_resolution(40, 4).band_rad_s, 0.01)
    };
    let smallest = (1..=40).find(|p| law(*p).is_ok()).expect("40 poles resolve 0.87 mm");
    assert!(smallest < 40, "the ceiling is not the need: {smallest}");
    let uniform = vec![Segment::Cylinder { radius, length: 0.57 }];
    let fixed = walled(uniform.clone(), 3, 4, 40, false);
    let adaptive = walled(uniform, 3, 4, 40, true);
    let memory = |d: &DgDuct| d.spatial_viscothermal_admission().unwrap().memory_dofs;
    assert_eq!(memory(&fixed), 2 * 3 * 4 * 40);
    assert_eq!(memory(&adaptive), 2 * 3 * 4 * smallest);
    // a ceiling below the need still refuses, with the same error as the fixed count
    let narrow = DgDuct::new(&Duct { segments: vec![Segment::Cylinder { radius: 0.019, length: 0.57 }] }, &g, DgMesh { elements: 3, order: 1 }, DgBoundary::Port(PortId(0)), DgBoundary::Closed).unwrap();
    assert!(narrow.with_spatial_viscothermal_min_poles(wall_resolution(8, 4)).is_err());
}

#[test]
fn min_poles_on_a_cone_lowers_the_dofs_and_the_stiffness_and_matches_the_dense_propagator() {
    // 4 elements, order 1, 2 wall points, ceiling 40 on a 0.87 mm to 19 mm cone: the thin end needs few poles.
    let fixed = walled(cone(), 4, 2, 40, false);
    let adaptive = walled(cone(), 4, 2, 40, true);
    let memory = |d: &DgDuct| d.spatial_viscothermal_admission().unwrap().memory_dofs;
    assert!(memory(&adaptive) * 10 < memory(&fixed) * 8, "{} against {}", memory(&adaptive), memory(&fixed));
    let (pf, pa) = (fixed.build_krylov_propagator(DT, TOL).unwrap(), adaptive.build_krylov_propagator(DT, TOL).unwrap());
    println!("KRYLOV-MEASURE cone E4: fixed {} DOFs stiffness {:.1}; min poles {} DOFs stiffness {:.1}", fixed.total_dofs(), stiffness(&pf), adaptive.total_dofs(), stiffness(&pa));
    assert!(stiffness(&pa) * 5.0 < stiffness(&pf));
    // the adaptive operator is still exactly what the dense propagator exponentiates
    assert_matches_dense(&adaptive, 200);
}

#[test]
fn min_poles_and_the_fixed_ceiling_agree_in_the_time_domain() {
    // Both resolve the same physics to the declared wall error, so the port mean pressure of a driven duct agrees to well
    // under that error. Compared on the Krylov kernel for both (the fixed duct is above what a debug dense build can do).
    let (fixed, adaptive) = (walled(cone(), 4, 2, 40, false), walled(cone(), 4, 2, 40, true));
    let (pf, pa) = (fixed.build_krylov_propagator(DT, TOL).unwrap(), adaptive.build_krylov_propagator(DT, TOL).unwrap());
    let (mut sf, mut sa) = (fixed.zero_state(), adaptive.zero_state());
    let (mut num, mut den) = (0.0_f64, 0.0_f64);
    for k in 0..200 {
        let a = pf.step(&mut sf, &input(effort(k))).unwrap().port_outputs[0].conjugate_rate;
        let b = pa.step(&mut sa, &input(effort(k))).unwrap().port_outputs[0].conjugate_rate;
        num += (a - b) * (a - b);
        den += a * a;
    }
    let difference = det::sqrt(num / den);
    println!("KRYLOV-MEASURE min poles against the fixed ceiling, port mean pressure over 200 steps: relative rms difference {difference:e}");
    assert!(difference < 1e-2, "{difference}");
}

#[test]
fn an_oboe_sized_bore_with_min_poles_is_small_stiffness_light_and_one_substep() {
    // The oboe bore of the test above, 39 elements, order 1, 2 wall points, with the ceiling 40 the bell needs.
    let duct = walled(oboe_like(), 39, 2, 40, true);
    assert!(duct.total_dofs() <= 4800, "{} DOFs", duct.total_dofs());
    let propagator = duct.propagator(DT, TOL).unwrap();
    let mut state = duct.zero_state();
    let mut previous = f64::INFINITY;
    for k in 0..40 {
        let response = propagator.port_response(&state, PortId(0)).unwrap();
        let u = effort(k);
        let out = propagator.step(&mut state, &input(u)).unwrap();
        let promised = response.history_pressure_pa + response.impedance_pa_s_m3 * u;
        let delivered = out.port_outputs[0].conjugate_rate;
        assert!((promised - delivered).abs() <= 1e-9 * (promised.abs() + delivered.abs()).max(1e-30), "step {k}: {promised} vs {delivered}");
        if k >= 24 {
            assert!(out.ledger.stored_energy_j <= previous * (1.0 + 1e-12) + 1e-300, "step {k}");
        }
        previous = out.ledger.stored_energy_j;
    }
    let report = propagator.krylov_report().expect("above the dense cap");
    println!("KRYLOV-MEASURE oboe-sized bore, min poles: {report:?}, stiffness {:.1}", stiffness(&propagator));
    assert_eq!(report.substeps, 1, "{report:?}");
    assert!(report.max_dimension <= 30, "{report:?}");
    assert!(stiffness(&propagator) <= 60.0, "{}", stiffness(&propagator));
}
