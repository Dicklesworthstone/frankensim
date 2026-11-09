//! G1/G3 analytic axial elasticity, compatibility, energy and load reversal.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_truss::elastic::{ElasticControls, ElasticError, solve};
use fs_truss::{GroundLimits, GroundStructure, LayoutCase};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey { seed: 1, kernel_id: 1, tile: 0, iteration: 0 },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-11 * a.abs().max(b.abs()).max(1e-15),
        "{a:e} != {b:e}"
    );
}

fn star(cx: &Cx<'_>, angle: f64, force: [f64; 2]) -> (GroundStructure, LayoutCase) {
    let (s, c) = angle.sin_cos();
    let rotate = |p: [f64; 2]| [c * p[0] - s * p[1], s * p[0] + c * p[1]];
    let nodes = [[-1.0, 0.0], [0.0, -1.0], [-1.0, -1.0], [0.0, 0.0]].map(rotate);
    let ground = GroundStructure::try_from_parts(
        &nodes,
        &[(0, 3), (1, 3), (2, 3)],
        &[1.0, 1.0, 2.0_f64.sqrt()],
        GroundLimits::default(),
        cx,
    )
    .unwrap();
    let case = LayoutCase::try_new(
        vec![[true; 2], [true; 2], [true; 2], [false; 2]],
        vec![[0.0; 2], [0.0; 2], [0.0; 2], rotate(force)],
        4,
    )
    .unwrap();
    (ground, case)
}

#[test]
fn g1_series_bars_keep_interface_compatibility_and_supported_loads() {
    with_cx(|cx, _| {
        let ground = GroundStructure::try_from_parts(
            &[[0.0, 0.0], [2.0, 0.0], [3.0, 0.0]],
            &[(0, 1), (1, 2)],
            &[2.0, 1.0],
            GroundLimits::default(),
            cx,
        )
        .unwrap();
        for sign in [-1.0, 1.0] {
            let case = LayoutCase::try_new(
                vec![[true, true], [false, true], [false, true]],
                vec![[2.0 * sign, -3.0], [0.0; 2], [12.0 * sign, 0.0]],
                3,
            )
            .unwrap();
            let response = solve(
                &ground, &case, &[2.0, 3.0], &[4e6, 2e6],
                ElasticControls::default(), cx,
            )
            .unwrap();
            near(response.displacements_m[1][0], 3e-6 * sign);
            near(response.displacements_m[2][0], 5e-6 * sign);
            near(response.member_force_n[0], 12.0 * sign);
            near(response.member_force_n[1], 12.0 * sign);
            near(response.member_stress_pa[0], 6.0 * sign);
            near(response.member_stress_pa[1], 4.0 * sign);
            near(response.support_reactions_n[0][0], -14.0 * sign);
            near(response.support_reactions_n[0][1], 3.0);
            near(response.strain_energy_j, 30e-6);
            near(response.compliance_j, 60e-6);
            assert!(response.relative_equilibrium_residual < 1e-12);
        }
    });
}

#[test]
fn g1_redundant_truss_shares_load_by_elastic_compatibility_not_equilibrium_alone() {
    with_cx(|cx, _| {
        let (ground, case) = star(cx, 0.0, [11.0, 15.0]);
        // Independently assembled K = 1e6 [[5,3],[3,6]], so u=1e-6 [1,2].
        let moduli = [2e6, 3e6, 6e6 * 2.0_f64.sqrt()];
        let response = solve(
            &ground, &case, &[1.0; 3], &moduli, ElasticControls::default(), cx,
        )
        .unwrap();
        near(response.displacements_m[3][0], 1e-6);
        near(response.displacements_m[3][1], 2e-6);
        near(response.member_force_n[0], 2.0);
        near(response.member_force_n[1], 6.0);
        near(response.member_force_n[2], 9.0 * 2.0_f64.sqrt());
        near(response.strain_energy_j, 20.5e-6);
        near(response.compliance_j, 41e-6);
        for axis in 0..2 {
            near(
                response.support_reactions_n.iter().map(|r| r[axis]).sum(),
                -case.loads()[3][axis],
            );
        }
        // [11,15,0] also balances loads, but cannot meet at one elastic node.
        assert!((response.member_force_n[0] - 11.0).abs() > 1.0);
        let stiffer = solve(
            &ground, &case, &[1.0, 1.0, 2.0], &moduli, ElasticControls::default(), cx,
        )
        .unwrap();
        near(stiffer.displacements_m[3][0], 0.25e-6);
        near(stiffer.displacements_m[3][1], 1.5e-6);
        near(stiffer.member_force_n[0], 0.5);
        near(stiffer.member_force_n[1], 4.5);
        near(stiffer.member_force_n[2], 10.5 * 2.0_f64.sqrt());
    });
}

#[test]
fn g3_rotating_the_stable_truss_preserves_stress_energy_and_rotates_motion() {
    with_cx(|cx, _| {
        for angle in [0.0, 0.37, 1.1] {
            let (ground, case) = star(cx, angle, [11.0, 15.0]);
            let response = solve(
                &ground, &case, &[1.0; 3], &[2e6, 3e6, 6e6 * 2.0_f64.sqrt()],
                ElasticControls::default(), cx,
            )
            .unwrap();
            let (s, c) = angle.sin_cos();
            near(response.displacements_m[3][0], 1e-6 * (c - 2.0 * s));
            near(response.displacements_m[3][1], 1e-6 * (s + 2.0 * c));
            near(response.member_force_n[0], 2.0);
            near(response.member_force_n[1], 6.0);
            near(response.strain_energy_j, 20.5e-6);
        }
    });
}

#[test]
fn g0_pruned_members_and_unloaded_coordinates_receive_no_artificial_stiffness() {
    with_cx(|cx, _| {
        let (ground, case) = star(cx, 0.0, [11.0, 15.0]);
        let response = solve(
            &ground, &case, &[1.0, 1.0, 0.0], &[2e6, 3e6, 1e9],
            ElasticControls::default(), cx,
        )
        .unwrap();
        near(response.displacements_m[3][0], 5.5e-6);
        near(response.displacements_m[3][1], 5e-6);
        assert_eq!(response.member_force_n[2], 0.0);
        assert_eq!(response.member_stress_pa[2], 0.0);
        let zero = LayoutCase::try_new(vec![[false; 2]; 4], vec![[0.0; 2]; 4], 4).unwrap();
        let response = solve(
            &ground, &zero, &[0.0; 3], &[1.0; 3], ElasticControls::default(), cx,
        )
        .unwrap();
        assert_eq!(response.displacements_m, vec![[0.0; 2]; 4]);
        assert_eq!(response.strain_energy_j, 0.0);
        assert_eq!(response.iterations, 0);
    });
}

#[test]
fn g4_unresolved_loads_invalid_materials_budget_and_cancellation_refuse() {
    with_cx(|cx, gate| {
        let (ground, case) = star(cx, 0.0, [11.0, 15.0]);
        let controls = ElasticControls::default();
        assert!(matches!(
            solve(&ground, &case, &[1.0, 0.0, 0.0], &[1.0; 3], controls, cx),
            Err(ElasticError::UnsupportedLoad { node: 3, axis: 1 })
        ));
        // A lone diagonal bar has a loaded transverse mechanism although
        // both Cartesian diagonal stiffness entries are positive.
        assert!(solve(&ground, &case, &[0.0, 0.0, 1.0], &[1.0; 3], controls, cx).is_err());
        assert!(solve(
            &ground, &case, &[1.0; 3], &[1.0; 3],
            ElasticControls { maximum_iterations: 0, ..controls }, cx,
        ).is_err());
        assert!(solve(&ground, &case, &[1.0; 2], &[1.0; 3], controls, cx).is_err());
        assert!(solve(&ground, &case, &[1.0; 3], &[1.0, f64::NAN, 1.0], controls, cx).is_err());
        gate.request();
        assert!(matches!(
            solve(&ground, &case, &[1.0; 3], &[1.0; 3], controls, cx),
            Err(ElasticError::Construction(fs_truss::TrussConstructionError::Cancelled { .. }))
        ));
    });
}
