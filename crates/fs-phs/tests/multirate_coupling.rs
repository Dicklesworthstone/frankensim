use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_phs::multirate::{
    Clock, IntegrationMethod, MultirateError, MultiratePair, PortPair, PowerTransfer, SolveConfig,
};
use fs_phs::{PortHamiltonian, QuadraticStorage, Storage, duffing_oscillator, mass_spring_damper};

fn exact(substeps: usize) -> Clock {
    Clock { substeps, method: IntegrationMethod::ExactQuadratic, maximum_contact_step_s: None }
}
fn nonlinear(substeps: usize) -> Clock {
    Clock { substeps, method: IntegrationMethod::DiscreteGradient {
        maximum_rate_rad_s: 5.0, maximum_rate_dt: 0.5 }, maximum_contact_step_s: None }
}
fn pair() -> [PortPair;1] { [PortPair { left:0, right:0, input_scale:1.0 }] }
fn close(a:f64,b:f64,tol:f64) { assert!((a-b).abs()<=tol,"{a:e} != {b:e}"); }
fn distance(a:&[f64],b:&[f64])->f64 {
    a.iter().zip(b).map(|(x,y)|(x-y)*(x-y)).sum::<f64>().sqrt()
}

#[test]
fn lossless_energy_remains_closed_for_ten_thousand_windows_on_different_clocks() {
    let left=mass_spring_damper(1.0,1.0,0.0).unwrap();
    let right=mass_spring_damper(1.0,1.0,0.0).unwrap();
    for substeps in [1,2,4,8] {
        let mut coupling=MultiratePair::new(&left,&right,&pair(),0.005,
            exact(1),exact(substeps),SolveConfig::default()).unwrap();
        let (mut a,mut b)=(vec![0.7,0.2],vec![-0.3,0.4]);
        let initial=left.hamiltonian(&a)+right.hamiltonian(&b);
        for _ in 0..10_000 {
            let r=coupling.step(&a,&b,&[0.0],&[0.0]).unwrap();
            close(r.exchange_work_j.iter().sum(),0.0,2e-13);
            close(r.balance_residual_j,0.0,2e-13);
            assert_eq!(r.external_work_j,[0.0,0.0]);
            assert_eq!(r.left.substeps,1);
            assert_eq!(r.right.substeps,substeps);
            a=r.left.state; b=r.right.state;
        }
        close(left.hamiltonian(&a)+right.hamiltonian(&b),initial,2e-10);
    }
}

#[test]
fn source_work_and_physical_damping_are_separate_from_interface_exchange() {
    let left=mass_spring_damper(1.0,2.0,0.3).unwrap();
    let right=mass_spring_damper(2.0,3.0,0.4).unwrap();
    let mut coupling=MultiratePair::new(&left,&right,&pair(),0.01,
        exact(2),exact(8),SolveConfig::default()).unwrap();
    let (mut a,mut b)=(vec![0.5,0.3],vec![0.1,-0.2]);
    let initial=left.hamiltonian(&a)+right.hamiltonian(&b);
    let (mut work,mut heat)=(0.0,0.0);
    for _ in 0..300 {
        let r=coupling.step(&a,&b,&[0.7],&[-0.2]).unwrap();
        assert!(r.left.dissipated_j>=0.0 && r.right.dissipated_j>=0.0);
        close(r.exchange_work_j.iter().sum(),0.0,1e-13);
        work+=r.external_work_j.iter().sum::<f64>();
        heat+=r.left.dissipated_j+r.right.dissipated_j;
        a=r.left.state; b=r.right.state;
    }
    assert!(heat>0.01);
    close(left.hamiltonian(&a)+right.hamiltonian(&b)-initial+heat-work,0.0,2e-11);
}

fn free_mass(mass: f64) -> PortHamiltonian {
    // Momentum is the state; force is input and physical velocity is output.
    PortHamiltonian::new(1,1,vec![0.0],vec![0.0],vec![1.0],
        Box::new(QuadraticStorage::new(vec![1.0/mass],1).unwrap())).unwrap()
}

fn inertias(masses: &[f64]) -> PortHamiltonian {
    let n = masses.len();
    let mut storage = vec![0.0; n * n];
    let mut ports = vec![0.0; n * n];
    for i in 0..n {
        storage[i * n + i] = 1.0 / masses[i];
        ports[i * n + i] = 1.0;
    }
    PortHamiltonian::new(
        n,
        n,
        vec![0.0; n * n],
        vec![0.0; n * n],
        ports,
        Box::new(QuadraticStorage::new(storage, n).unwrap()),
    )
    .unwrap()
}

fn piston_transfer() -> PowerTransfer {
    PowerTransfer {
        left_ports: vec![0, 1],
        right_ports: vec![0],
        right_from_left: vec![0.1, -0.07],
        left_input_scales: vec![1.0, 1.0],
    }
}

#[test]
fn unequal_signed_piston_areas_drive_one_compliance_with_reciprocal_force_and_work() {
    // G0/G1: two inertias exchange power with ONE compressible chamber.
    // p_i' = F_i - A_i P, V' = Q + sum A_i p_i/m_i, P = V/C.
    // The opposite piston orientation is a signed physical area, not a
    // numerical scale. Solve these window-average equations independently.
    let masses = [2.0, 3.0];
    let compliance = 0.004;
    let (left, right) = (inertias(&masses), free_mass(compliance));
    let transfer = piston_transfer();
    let (p, v, f, q, dt) = ([0.8, -0.3], 0.003, [0.3, -0.2], 0.01, 0.2);
    let areas = [0.1, -0.07];
    let stiffness: f64 = (0..2).map(|i| areas[i] * areas[i] / masses[i]).sum();
    let free_flow: f64 = (0..2).map(|i| areas[i] * p[i] / masses[i]).sum();
    let forced_flow: f64 = (0..2).map(|i| areas[i] * f[i] / masses[i]).sum();
    let pressure = (v / compliance
        + dt / (2.0 * compliance) * (q + free_flow)
        + dt * dt / (4.0 * compliance) * forced_flow)
        / (1.0 + dt * dt / (4.0 * compliance) * stiffness);
    let momentum: Vec<f64> = (0..2)
        .map(|i| p[i] + dt * (f[i] - areas[i] * pressure))
        .collect();
    let velocity: Vec<f64> = (0..2)
        .map(|i| (p[i] + momentum[i]) / (2.0 * masses[i]))
        .collect();
    let flow: f64 = (0..2).map(|i| areas[i] * velocity[i]).sum();
    let mut coupling = MultiratePair::new_mapped(
        &left,
        &right,
        &transfer,
        dt,
        exact(2),
        exact(8),
        SolveConfig::default(),
    )
    .unwrap();
    let r = coupling.step(&p, &[v], &f, &[q]).unwrap();
    for i in 0..2 {
        close(r.left.state[i], momentum[i], 2e-14);
        close(r.left.average_output[i], velocity[i], 2e-14);
    }
    close(r.right.state[0], v + dt * (q + flow), 2e-14);
    close(r.right.average_output[0], pressure, 2e-13);
    close(
        r.external_work_j[0],
        dt * (f[0] * velocity[0] + f[1] * velocity[1]),
        2e-14,
    );
    close(r.external_work_j[1], dt * q * pressure, 2e-14);
    close(r.exchange_work_j[0], -dt * pressure * flow, 2e-14);
    close(r.exchange_work_j[1], dt * pressure * flow, 2e-14);
    close(r.balance_residual_j, 0.0, 2e-14);
    assert!(r.exchange_work_j[0] < 0.0 && r.exchange_work_j[1] > 0.0);
    // The same physical assembly with the two owners exchanged uses -M^T.
    // This also exercises a tall (two-by-one) map and independent solve scaling.
    let reversed = PowerTransfer {
        left_ports: vec![0],
        right_ports: vec![0, 1],
        right_from_left: vec![-0.1, 0.07],
        left_input_scales: vec![3.0],
    };
    let s = MultiratePair::new_mapped(
        &right,
        &left,
        &reversed,
        dt,
        exact(8),
        exact(2),
        SolveConfig::default(),
    )
    .unwrap()
    .step(&[v], &p, &[q], &f)
    .unwrap();
    assert!(distance(&r.left.state, &s.right.state) < 2e-14);
    assert!(distance(&r.right.state, &s.left.state) < 2e-14);
    for i in 0..2 {
        close(r.external_work_j[i], s.external_work_j[1 - i], 2e-14);
        close(r.exchange_work_j[i], s.exchange_work_j[1 - i], 2e-14);
    }
}

#[test]
fn rectangular_transfer_retains_null_motion_and_converges_to_physical_coupled_frequency() {
    // G1/G3: the piston system has one oscillating pressure mode and one
    // uncoupled momentum combination. Neither may be removed by an SVD cutoff.
    let masses = [2.0, 3.0];
    let compliance = 0.004;
    let areas = [0.1, -0.07];
    let (left, right) = (inertias(&masses), free_mass(compliance));
    let transfer = piston_transfer();
    let initial = [0.8, -0.3];
    let volume = 0.003;
    let duration = 2.0;
    let omega: f64 =
        ((areas[0] * areas[0] / masses[0] + areas[1] * areas[1] / masses[1]) / compliance).sqrt();
    let flow = areas[0] * initial[0] / masses[0] + areas[1] * initial[1] / masses[1];
    let angle = omega * duration;
    let volume_integral =
        volume * angle.sin() / omega + flow * (1.0 - angle.cos()) / (omega * omega);
    let truth = [
        initial[0] - areas[0] * volume_integral / compliance,
        initial[1] - areas[1] * volume_integral / compliance,
        volume * angle.cos() + flow * angle.sin() / omega,
    ];
    let mut errors = Vec::new();
    for windows in [10, 20, 40] {
        let mut coupling = MultiratePair::new_mapped(
            &left,
            &right,
            &transfer,
            duration / windows as f64,
            exact(2),
            exact(8),
            SolveConfig::default(),
        )
        .unwrap();
        let (mut p, mut v) = (initial.to_vec(), vec![volume]);
        let energy = left.hamiltonian(&p) + right.hamiltonian(&v);
        for _ in 0..windows {
            let r = coupling.step(&p, &v, &[0.0; 2], &[0.0]).unwrap();
            close(r.exchange_work_j.iter().sum(), 0.0, 2e-14);
            p = r.left.state;
            v = r.right.state;
        }
        close(
            areas[1] * p[0] - areas[0] * p[1],
            areas[1] * initial[0] - areas[0] * initial[1],
            2e-14,
        );
        close(left.hamiltonian(&p) + right.hamiltonian(&v), energy, 2e-13);
        p.extend(v);
        errors.push(distance(&p, &truth));
    }
    assert!(
        errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8,
        "{errors:?}"
    );
}

#[test]
fn mapped_interfaces_preserve_port_ordering_scales_and_unconnected_sources() {
    // G3: a full non-diagonal map and its reordered coordinates describe the
    // same physics. The third native port stays externally driven on each side.
    let (left, right) = (two_modes(), two_modes());
    let map = PowerTransfer {
        left_ports: vec![0, 1],
        right_ports: vec![0, 1],
        right_from_left: vec![0.5, -2.0, 1.2, 0.3],
        left_input_scales: vec![1.0, 1.0],
    };
    let reordered = PowerTransfer {
        left_ports: vec![1, 0],
        right_ports: vec![1, 0],
        right_from_left: vec![0.3, 1.2, -2.0, 0.5],
        left_input_scales: vec![100.0, 0.01],
    };
    let solve = |map: &PowerTransfer| {
        MultiratePair::new_mapped(
            &left,
            &right,
            map,
            0.03,
            exact(2),
            exact(8),
            SolveConfig::default(),
        )
        .unwrap()
        .step_subcycled_sources(
            &[0.3, 0.4, -0.2, 0.1],
            &[-0.1, 0.2, 0.7, -0.3],
            &[0.4, -0.1, 0.2, -0.2, 0.3, -0.1],
            &[-0.2, 0.3, -0.1].repeat(8),
        )
        .unwrap()
    };
    let (a, b) = (solve(&map), solve(&reordered));
    assert!(distance(&a.left.state, &b.left.state) < 1e-13);
    assert!(distance(&a.right.state, &b.right.state) < 1e-13);
    for i in 0..2 {
        close(a.external_work_j[i], b.external_work_j[i], 2e-14);
        close(a.exchange_work_j[i], b.exchange_work_j[i], 2e-14);
    }
    close(a.exchange_work_j.iter().sum(), 0.0, 2e-14);
}

#[test]
fn mapped_nonlinear_ports_keep_original_heat_and_cancel_without_changing_retry() {
    // G0/G4: the physical transfer acts on original nonlinear outputs; it does
    // not replace the material storage or loss law with a quadratic surrogate.
    let left = duffing_oscillator(1.0, 1.0, 0.8, 0.05).unwrap();
    let right = duffing_oscillator(1.2, 1.3, 0.6, 0.02).unwrap();
    let map = PowerTransfer {
        left_ports: vec![0],
        right_ports: vec![0],
        right_from_left: vec![-2.5],
        left_input_scales: vec![0.1],
    };
    let make = || {
        MultiratePair::new_mapped(
            &left,
            &right,
            &map,
            0.02,
            nonlinear(2),
            nonlinear(4),
            SolveConfig::default(),
        )
        .unwrap()
    };
    let (a, b) = ([0.6, 0.2], [-0.4, 0.1]);
    let mut candidate = make();
    let mut calls = 0;
    assert!(matches!(
        candidate.step_controlled(&a, &b, &[0.1], &[-0.05], || {
            calls += 1;
            calls >= 9
        }),
        Err(MultirateError::Cancelled)
    ));
    let r = candidate.step(&a, &b, &[0.1], &[-0.05]).unwrap();
    assert_eq!(r, make().step(&a, &b, &[0.1], &[-0.05]).unwrap());
    assert!(r.left.dissipated_j > 0.0 && r.right.dissipated_j > 0.0);
    close(r.exchange_work_j.iter().sum(), 0.0, 2e-13);
    close(r.balance_residual_j, 0.0, 2e-13);
}

#[test]
fn zero_map_retains_independent_evolution_and_invalid_maps_refuse_before_preparation() {
    let (left, right) = (inertias(&[2.0, 3.0]), free_mass(0.004));
    let zero = PowerTransfer {
        right_from_left: vec![0.0, 0.0],
        ..piston_transfer()
    };
    let r = MultiratePair::new_mapped(
        &left,
        &right,
        &zero,
        0.2,
        exact(2),
        exact(8),
        SolveConfig::default(),
    )
    .unwrap()
    .step(&[0.8, -0.3], &[0.003], &[0.3, -0.2], &[0.01])
    .unwrap();
    close(r.left.state[0], 0.86, 2e-14);
    close(r.left.state[1], -0.34, 2e-14);
    close(r.right.state[0], 0.005, 2e-14);
    assert_eq!(r.exchange_work_j, [0.0, 0.0]);
    let valid = piston_transfer();
    for bad in [
        PowerTransfer {
            left_ports: vec![],
            ..valid.clone()
        },
        PowerTransfer {
            left_ports: vec![0, 0],
            ..valid.clone()
        },
        PowerTransfer {
            right_ports: vec![1],
            ..valid.clone()
        },
        PowerTransfer {
            right_from_left: vec![0.1],
            ..valid.clone()
        },
        PowerTransfer {
            right_from_left: vec![0.1, f64::NAN],
            ..valid.clone()
        },
        PowerTransfer {
            left_input_scales: vec![1.0, 0.0],
            ..valid.clone()
        },
        PowerTransfer {
            left_input_scales: vec![1.0],
            ..valid.clone()
        },
        PowerTransfer {
            left_ports: vec![0; 33],
            ..valid
        },
    ] {
        assert!(matches!(
            MultiratePair::new_mapped(
                &left,
                &right,
                &bad,
                0.2,
                exact(2),
                exact(8),
                SolveConfig::default()
            ),
            Err(MultirateError::Invalid(_))
        ));
    }
}

#[test]
fn switched_native_forces_keep_their_work_even_when_both_mean_forces_vanish() {
    // G1: two inertias, exact constant-acceleration segments. Solving the
    // held interface independently gives cL=0.36 N, cR=1.59 N and mean
    // velocities 1.59 and -0.36 m/s. Source work is sum F_k * delta_x_k;
    // both average sources are zero, but their actual work is not.
    let (left,right)=(free_mass(2.0),free_mass(3.0));
    let mut coupling=MultiratePair::new(&left,&right,&pair(),1.0,
        exact(2),exact(4),SolveConfig::default()).unwrap();
    let r=coupling.step_subcycled_sources(&[2.0],&[-1.5],
        &[4.0,-4.0],&[-3.0,0.0,3.0,0.0]).unwrap();
    close(r.left.state[0],2.36,3e-14);
    close(r.right.state[0],0.09,3e-14);
    close(r.left.average_output[0],1.59,3e-14);
    close(r.right.average_output[0],-0.36,3e-14);
    close(r.external_work_j[0],-0.18,3e-14);
    close(r.external_work_j[1],0.19875,3e-14);
    close(r.exchange_work_j[0],0.5724,3e-14);
    close(r.exchange_work_j[1],-0.5724,3e-14);
    close(r.left.stored_change_j+r.right.stored_change_j,0.01875,3e-14);
    close(r.balance_residual_j,0.0,3e-14);
    assert_eq!(r.left.dissipated_j,0.0);
    assert_eq!(r.right.dissipated_j,0.0);
}

#[test]
fn repeated_native_sources_match_held_sources_bitwise_with_multiple_ports() {
    // G3: repeating a held input does not change its meaning, state or work.
    let (left,right)=(two_modes(),two_modes());
    let ports=[PortPair{left:0,right:1,input_scale:0.1},
        PortPair{left:1,right:0,input_scale:10.0}];
    let make=||MultiratePair::new(&left,&right,&ports,0.03,
        exact(2),exact(8),SolveConfig::default()).unwrap();
    let (a,b)=([0.3,0.4,-0.2,0.1],[-0.1,0.2,0.7,-0.3]);
    let (ul,ur)=([0.4,-0.1,0.2],[-0.2,0.3,-0.1]);
    let held=make().step(&a,&b,&ul,&ur).unwrap();
    let scheduled=make().step_subcycled_sources(&a,&b,&ul.repeat(2),&ur.repeat(8)).unwrap();
    assert_eq!(held,scheduled);
}

#[test]
fn nonlinear_subcycled_sources_close_physical_heat_and_regenerative_work() {
    // G0: varying source histories use the original nonlinear storage and
    // resistance at every native step; no work is fabricated from the balance.
    let left=duffing_oscillator(1.0,1.0,0.8,0.05).unwrap();
    let right=duffing_oscillator(1.2,1.3,0.6,0.02).unwrap();
    let mut coupling=MultiratePair::new(&left,&right,&pair(),0.02,
        nonlinear(2),nonlinear(4),SolveConfig::default()).unwrap();
    let (mut a,mut b)=(vec![0.6,0.2],vec![-0.4,0.1]);
    let initial=left.hamiltonian(&a)+right.hamiltonian(&b);
    let (mut work,mut heat)=(0.0,0.0);
    for _ in 0..50 {
        let r=coupling.step_subcycled_sources(&a,&b,&[0.1,-0.2],
            &[-0.05,0.2,0.0,0.1]).unwrap();
        close(r.exchange_work_j.iter().sum(),0.0,2e-13);
        close(r.balance_residual_j,0.0,2e-13);
        assert!(r.left.dissipated_j>=0.0 && r.right.dissipated_j>=0.0);
        work+=r.external_work_j.iter().sum::<f64>();
        heat+=r.left.dissipated_j+r.right.dissipated_j;
        a=r.left.state; b=r.right.state;
    }
    assert!(heat>0.001);
    close(left.hamiltonian(&a)+right.hamiltonian(&b)-initial+heat-work,0.0,2e-11);
}

#[test]
fn invalid_and_cancelled_source_histories_leave_the_interface_guess_reusable() {
    // G4: even cancellation after native substeps cannot publish partial state
    // or contaminate the next interface solve with a tentative source history.
    let (left,right)=(free_mass(2.0),free_mass(3.0));
    let make=||MultiratePair::new(&left,&right,&pair(),1.0,
        exact(2),exact(4),SolveConfig::default()).unwrap();
    let mut coupling=make();
    let (a,b)=([2.0],[-1.5]);
    let (ul,ur)=([4.0,-4.0],[-3.0,0.0,3.0,0.0]);
    for bad in [&[4.0][..],&[4.0,-4.0,0.0][..],&[4.0,f64::NAN][..]] {
        assert!(matches!(coupling.step_subcycled_sources(&a,&b,bad,&ur),
            Err(MultirateError::Invalid(_))));
    }
    assert!(matches!(coupling.step_subcycled_sources(&a,&b,&ul,&ur[..3]),
        Err(MultirateError::Invalid(_))));
    let mut calls=0;
    assert!(matches!(coupling.step_subcycled_sources_controlled(&a,&b,&ul,&ur,||{
        calls+=1;calls>=7
    }),Err(MultirateError::Cancelled)));
    assert_eq!(coupling.step_subcycled_sources(&a,&b,&ul,&ur).unwrap(),
        make().step_subcycled_sources(&a,&b,&ul,&ur).unwrap());
}

fn monolithic() -> PortHamiltonian {
    // Original states [qL,pL,qR,pR], including the full reciprocal reaction.
    PortHamiltonian::new(4,0,
        vec![0.0,1.0,0.0,0.0, -1.0,0.0,0.0,-1.0,
            0.0,0.0,0.0,1.0, 0.0,1.0,-1.0,0.0],vec![0.0;16],vec![],
        Box::new(QuadraticStorage::new(vec![1.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0,
            0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0],4).unwrap())).unwrap()
}

#[test]
fn refining_the_common_window_converges_to_the_full_coupled_exact_flow() {
    let left=mass_spring_damper(1.0,1.0,0.0).unwrap();
    let right=mass_spring_damper(1.0,1.0,0.0).unwrap();
    let whole=monolithic();
    let truth=PreparedExactZoh::new(&whole,2.0,LinearSplitBudget {
        max_modes:4,max_omega_dt:0.04 }).unwrap().step(&[0.7,0.2,-0.3,0.4],&[]).unwrap().x;
    let mut errors=Vec::new();
    for windows in [20,40,80] {
        let mut coupling=MultiratePair::new(&left,&right,&pair(),2.0/windows as f64,
            exact(1),exact(4),SolveConfig::default()).unwrap();
        let (mut a,mut b)=(vec![0.7,0.2],vec![-0.3,0.4]);
        for _ in 0..windows {
            let r=coupling.step(&a,&b,&[0.0],&[0.0]).unwrap();
            a=r.left.state; b=r.right.state;
        }
        a.extend(b); errors.push(distance(&a,&truth));
    }
    assert!(errors[0]/errors[1]>3.8,"{errors:?}");
    assert!(errors[1]/errors[2]>3.8,"{errors:?}");
}

#[test]
fn nonlinear_native_subcycles_resolve_the_original_duffing_law() {
    let left=duffing_oscillator(1.0,1.0,0.8,0.05).unwrap();
    let right=duffing_oscillator(1.2,1.3,0.6,0.02).unwrap();
    let evolve=|substeps| {
        let mut coupling=MultiratePair::new(&left,&right,&pair(),0.04,
            nonlinear(substeps),nonlinear(substeps),SolveConfig::default()).unwrap();
        let (mut a,mut b)=(vec![0.6,0.2],vec![-0.4,0.1]);
        for _ in 0..25 {
            let r=coupling.step(&a,&b,&[0.1],&[-0.05]).unwrap();
            close(r.balance_residual_j,0.0,2e-12);
            a=r.left.state; b=r.right.state;
        }
        a.extend(b); a
    };
    let reference=evolve(16);
    let errors:Vec<f64>=[1,2,4,8].into_iter().map(|n|distance(&evolve(n),&reference)).collect();
    assert!(errors.windows(2).all(|e|e[0]>3.0*e[1]),"{errors:?}");
}

fn two_modes() -> PortHamiltonian {
    PortHamiltonian::new(4,3,
        vec![0.0,1.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,-1.0,0.0],
        vec![0.0;16],vec![0.0,0.0,0.3, 1.0,0.0,0.0, 0.0,0.0,0.0, 0.0,1.0,0.0],
        Box::new(QuadraticStorage::new(vec![1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,
            0.0,0.0,2.0,0.0,0.0,0.0,0.0,0.5],4).unwrap())).unwrap()
}

#[test]
fn multiple_ports_keep_independent_sources_and_invariant_numerical_scaling() {
    let (left,right)=(two_modes(),two_modes());
    let normal=[PortPair{left:0,right:1,input_scale:1.0},PortPair{left:1,right:0,input_scale:1.0}];
    let scaled=[PortPair{left:1,right:0,input_scale:100.0},PortPair{left:0,right:1,input_scale:0.01}];
    let solve=|pairs:&[PortPair]| {
        MultiratePair::new(&left,&right,pairs,0.03,exact(2),exact(8),SolveConfig::default())
            .unwrap().step(&[0.3,0.4,-0.2,0.1],&[-0.1,0.2,0.7,-0.3],
                &[0.4,-0.1,0.2],&[-0.2,0.3,-0.1]).unwrap()
    };
    let (a,b)=(solve(&normal),solve(&scaled));
    assert!(distance(&a.left.state,&b.left.state)<1e-13);
    assert!(distance(&a.right.state,&b.right.state)<1e-13);
    close(a.external_work_j.iter().sum(),b.external_work_j.iter().sum(),1e-13);
    assert!(a.external_work_j.iter().any(|x|x.abs()>1e-6));
    close(a.exchange_work_j.iter().sum(),0.0,1e-13);
}

#[test]
fn clock_contact_and_iteration_limits_are_explicit() {
    let (left,right)=(mass_spring_damper(1.0,1.0,0.0).unwrap(),mass_spring_damper(1.0,1.0,0.0).unwrap());
    for n in [0,17] {
        assert!(MultiratePair::new(&left,&right,&pair(),0.01,exact(n),exact(1),SolveConfig::default()).is_err());
    }
    let unresolved=Clock{maximum_contact_step_s:Some(0.001),..exact(1)};
    assert!(MultiratePair::new(&left,&right,&pair(),0.01,unresolved,exact(1),SolveConfig::default()).is_err());
    assert!(MultiratePair::new(&left,&right,&pair(),1.0,nonlinear(1),exact(1),SolveConfig::default()).is_err());
    let mut zero=MultiratePair::new(&left,&right,&pair(),0.01,exact(1),exact(1),SolveConfig {
        max_iterations:0,..SolveConfig::default() }).unwrap();
    assert!(matches!(zero.step(&[0.1,0.3],&[0.2,0.4],&[0.0],&[0.0]),Err(MultirateError::NotConverged{..})));
}

#[test]
fn cancellation_publishes_nothing_and_retry_matches_a_fresh_coupler() {
    let (left,right)=(mass_spring_damper(1.0,1.0,0.0).unwrap(),mass_spring_damper(1.0,1.0,0.0).unwrap());
    let make=||MultiratePair::new(&left,&right,&pair(),0.02,exact(2),exact(8),SolveConfig::default()).unwrap();
    let mut candidate=make();
    let mut calls=0;
    let failed=candidate.step_controlled(&[0.4,0.1],&[0.2,-0.2],&[0.3],&[0.1],||{
        calls+=1; calls>=7
    });
    assert!(matches!(failed,Err(MultirateError::Cancelled)));
    let retry=candidate.step(&[0.4,0.1],&[0.2,-0.2],&[0.3],&[0.1]).unwrap();
    let fresh=make().step(&[0.4,0.1],&[0.2,-0.2],&[0.3],&[0.1]).unwrap();
    assert_eq!(retry,fresh);
}

#[test]
fn prepared_stiff_loss_integral_matches_the_analytic_held_input_solution() {
    let system=PortHamiltonian::new(1,1,vec![0.0],vec![100.0],vec![1.0],
        Box::new(QuadraticStorage::new(vec![1.0],1).unwrap())).unwrap();
    let r=PreparedExactZoh::new(&system,0.1,LinearSplitBudget{max_modes:1,max_omega_dt:0.04})
        .unwrap().step(&[1.0],&[1.0]).unwrap();
    let decay=(-10.0_f64).exp();
    let endpoint=0.01+0.99*decay;
    let work=0.001+0.99*(1.0-decay)/100.0;
    let heat=100.0*(0.0001*0.1+2.0*0.01*0.99*(1.0-decay)/100.0
        +0.99*0.99*(1.0-decay*decay)/200.0);
    close(r.x[0],endpoint,2e-14); close(r.supplied,work,2e-14);
    close(r.dissipated,heat,2e-13); close(r.delta_h,0.5*(endpoint*endpoint-1.0),2e-13);
}

struct DatumStorage { datum:f64, defect:f64 }
impl Storage for DatumStorage {
    fn hamiltonian(&self,x:&[f64])->f64 { self.datum+0.5*x[0]*x[0] }
    fn gradient(&self,x:&[f64],out:&mut[f64]) { out[0]=x[0]; }
    fn midpoint_energy_remainder(&self,a:&[f64],b:&[f64],mid:&[f64])->Option<f64> {
        let dx=b[0]-a[0];
        Some(0.5*dx*((a[0]-mid[0])+(b[0]-mid[0]))+self.defect*dx*dx*dx)
    }
}

#[test]
fn energy_admission_cannot_be_diluted_by_an_arbitrary_hamiltonian_datum() {
    for datum in [0.0,1e30] {
        for defect in [0.0,100.0] {
            let law=||PortHamiltonian::new(1,1,vec![0.0],vec![0.0],vec![1.0],
                Box::new(DatumStorage{datum,defect})).unwrap();
            let (left,right)=(law(),law());
            let mut coupling=MultiratePair::new(&left,&right,&pair(),0.04,
                nonlinear(2),nonlinear(4),SolveConfig::default()).unwrap();
            let result=coupling.step(&[0.7],&[0.8],&[0.2],&[0.3]);
            if defect==0.0 { assert!(result.is_ok(),"{result:?}"); }
            else { assert!(matches!(result,Err(MultirateError::Energy{..})),"{result:?}"); }
        }
    }
}
