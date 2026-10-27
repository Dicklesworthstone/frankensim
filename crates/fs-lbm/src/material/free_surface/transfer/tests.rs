use super::*;

mod enthalpy;
use crate::d3q19::freesurface3::ContactModel3;
use crate::material::free_surface::plan_free_surface_from_material;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::surface::{
    LiquidCapillarity, ResolvedLiquidSurfaceState, resolve_liquid_surface_state,
};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::Dims;

fn with_cx<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 19,
                kernel_id: 28,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn liquid(rho: f64) -> ResolvedLiquidSurfaceState {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("density", Dims([-3, 1, 0, 0, 0, 0]), rho),
        ("dynamic_viscosity", Dims([-1, 1, -1, 0, 0, 0]), rho * 1e-4),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 500.0, 500.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic conservative transfer fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "test liquid".into(),
            phase: "liquid".into(),
            process: "analytic fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_liquid_surface_state(
        &card,
        &QueryPoint::new().with("T", 500.0).unwrap(),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn limits() -> TetLiquidTransferLimits {
    TetLiquidTransferLimits {
        max_tetrahedra: 8,
        max_overlap_tests: 10_000,
        max_wet_cells: 4096,
        relative_volume_tolerance: 1e-11,
        mass_tolerance_kg: 1e-14,
        momentum_tolerance_kg_m_s: 1e-16,
        angular_momentum_tolerance_kg_m2_s: 1e-18,
        maximum_centroid_shift_m: 1e-10,
        maximum_projection_energy_change_j: 1e-6,
    }
}
fn domain() -> MaterialFreeSurfaceDomain3<'static> {
    MaterialFreeSurfaceDomain3 {
        cells: [16; 3],
        gravity_m_s2: [0.0; 3],
        contact: ContactModel3::Neutral,
        solid_tiles: &[],
        max_cells: 4096,
    }
}
fn cube(origin: [f64; 3]) -> (Vec<[f64; 3]>, Vec<[u32; 4]>, Vec<[f64; 3]>) {
    let points: Vec<[f64; 3]> = (0..8)
        .map(|i| {
            core::array::from_fn(|a| {
                origin[a] + 0.002 + if i & (1 << a) == 0 { 0.0 } else { 0.004 }
            })
        })
        .collect();
    let velocities = points
        .iter()
        .map(|&p| {
            add(
                [0.01, -0.005, 0.002],
                cross([0.0, 0.0, 3.0], sub(p, add(origin, [0.004; 3]))),
            )
        })
        .collect();
    (
        points,
        vec![
            [0, 1, 3, 7],
            [0, 3, 2, 7],
            [0, 2, 6, 7],
            [0, 6, 4, 7],
            [0, 4, 5, 7],
            [0, 5, 1, 7],
        ],
        velocities,
    )
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "actual={a:e}, expected={b:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn spinning_translating_cube_preserves_analytic_moments_and_refines_projection_energy() {
    let source = liquid(1000.0);
    let (p, t, v) = cube([0.0; 3]);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mass = 1000.0 * 0.004_f64.powi(3);
        let translation = [0.01, -0.005, 0.002];
        let inertia = mass * 0.004_f64.powi(2) / 6.0;
        let analytic_p = scale(translation, mass);
        let analytic_l = add(cross([0.004; 3], analytic_p), [0.0, 0.0, inertia * 3.0]);
        let kinetic = 0.5 * mass * dot(translation, translation) + 0.5 * inertia * 9.0;
        let mut previous_error = f64::INFINITY;
        for dx in [0.001, 0.0005] {
            let scaling =
                plan_free_surface_from_material(&source, dx, dx * dx / 0.001, 0.1).unwrap();
            let (mut flow, report) = scaling
                .free_surface3_from_tetrahedra(
                    domain(),
                    [0.0; 3],
                    TetrahedralLiquidState {
                        positions_m: &p,
                        tetrahedra: &t,
                        velocities_m_s: &v,
                    },
                    limits(),
                    cx,
                )
                .unwrap();
            close(report.source.mass_kg, mass, 2e-18);
            close(report.target.mass_kg, mass, 2e-18);
            close(report.source.kinetic_energy_j, kinetic, 2e-22);
            for a in 0..3 {
                close(report.source.momentum_kg_m_s[a], analytic_p[a], 2e-20);
                close(report.target.momentum_kg_m_s[a], analytic_p[a], 2e-20);
                close(
                    report.source.angular_momentum_kg_m2_s[a],
                    analytic_l[a],
                    2e-22,
                );
                close(
                    report.target.angular_momentum_kg_m2_s[a],
                    analytic_l[a],
                    2e-22,
                );
                close(report.centroid_shift_m[a], 0.0, 2e-16);
            }
            // Point masses at regular voxel centres omit each cell's own
            // moment of inertia. Conserving spin therefore increases resolved
            // kinetic energy by this independently calculated projection error.
            let discrete_inertia = inertia * (1.0 - (dx / 0.004).powi(2));
            let expected_change = 0.5 * 9.0 * inertia * (inertia / discrete_inertia - 1.0);
            close(report.kinetic_energy_change_j, expected_change, 2e-22);
            assert!(report.kinetic_energy_change_j > 0.0);
            assert!(report.kinetic_energy_change_j < previous_error);
            previous_error = report.kinetic_energy_change_j;
            // The transferred field is a real native solver state, not a report.
            let initial_mass = flow.ledger_mass_kg().unwrap();
            flow.step().unwrap();
            close(flow.ledger_mass_kg().unwrap(), initial_mass, 1e-14);
        }
    });
}

#[test]
fn material_density_swap_scales_conserved_quantities_and_keeps_native_similitude() {
    let (p, t, v) = cube([0.0; 3]);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut previous = None;
        for rho in [1000.0, 2000.0] {
            let source = liquid(rho);
            let scaling = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1).unwrap();
            let (mut flow, report) = scaling
                .free_surface3_from_tetrahedra(
                    domain(),
                    [0.0; 3],
                    TetrahedralLiquidState {
                        positions_m: &p,
                        tetrahedra: &t,
                        velocities_m_s: &v,
                    },
                    limits(),
                    cx,
                )
                .unwrap();
            flow.step().unwrap();
            let result = (
                report.source.mass_kg,
                report.target.kinetic_energy_j,
                flow.lattice().grid().state_bits(),
            );
            if let Some((mass, energy, ref bits)) = previous {
                close(result.0, 2.0 * mass, 1e-18);
                close(result.1, 2.0 * energy, 1e-22);
                assert_eq!(&result.2, bits);
            }
            previous = Some(result);
        }
    });
}

#[test]
fn a_cut_tetrahedron_transfers_true_volume_and_affine_motion() {
    let source = liquid(1000.0);
    let p = [
        [0.0023, 0.0022, 0.0021],
        [0.0057, 0.0022, 0.0021],
        [0.0023, 0.0056, 0.0021],
        [0.0023, 0.0022, 0.0055],
    ];
    let v = p.map(|x| [0.005 + x[0], 2.0 * x[1], -x[2]]);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1).unwrap();
        let mut bounds = limits();
        bounds.maximum_centroid_shift_m = 0.001;
        let (_, report) = scaling
            .free_surface3_from_tetrahedra(
                domain(),
                [0.0; 3],
                TetrahedralLiquidState {
                    positions_m: &p,
                    tetrahedra: &[[0, 1, 2, 3]],
                    velocities_m_s: &v,
                },
                bounds,
                cx,
            )
            .unwrap();
        close(
            report.source.mass_kg,
            1000.0 * 0.0034_f64.powi(3) / 6.0,
            1e-18,
        );
        close(report.target.mass_kg, report.source.mass_kg, 1e-18);
        assert!(report.wet_cells > 1);
        assert!(report.centroid_shift_m.iter().any(|x| x.abs() > 1e-8));
        for a in 0..3 {
            close(
                report.source.momentum_kg_m_s[a],
                report.target.momentum_kg_m_s[a],
                1e-19,
            );
            close(
                report.source.angular_momentum_kg_m2_s[a],
                report.target.angular_momentum_kg_m2_s[a],
                1e-21,
            );
        }
    });
}

#[test]
fn insufficient_resolution_energy_or_work_and_cancellation_refuse_before_publication() {
    let source = liquid(1000.0);
    let (p, t, v) = cube([0.0; 3]);
    let gate = CancelGate::new_clock_free();
    let snapshot = TetrahedralLiquidState {
        positions_m: &p,
        tetrahedra: &t,
        velocities_m_s: &v,
    };
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1).unwrap();
        let mut no_energy_error = limits();
        no_energy_error.maximum_projection_energy_change_j = 0.0;
        assert!(
            scaling
                .free_surface3_from_tetrahedra(domain(), [0.0; 3], snapshot, no_energy_error, cx)
                .is_err()
        );
        let mut no_work = limits();
        no_work.max_overlap_tests = 1;
        assert!(
            scaling
                .free_surface3_from_tetrahedra(domain(), [0.0; 3], snapshot, no_work, cx)
                .is_err()
        );
        let mut no_cells = limits();
        no_cells.max_wet_cells = 1;
        assert!(
            scaling
                .free_surface3_from_tetrahedra(domain(), [0.0; 3], snapshot, no_cells, cx)
                .is_err()
        );
        // One target point cannot carry a tetrahedron's intrinsic angular
        // momentum. A generous energy/centroid budget must not hide this rank loss.
        let small_p = p
            .iter()
            .map(|x| x.map(|a| 0.002 + a * 0.01))
            .collect::<Vec<_>>();
        let mut coarse = limits();
        coarse.maximum_centroid_shift_m = 1.0;
        assert!(
            scaling
                .free_surface3_from_tetrahedra(
                    domain(),
                    [0.0; 3],
                    TetrahedralLiquidState {
                        positions_m: &small_p,
                        tetrahedra: &t,
                        velocities_m_s: &v
                    },
                    coarse,
                    cx
                )
                .is_err()
        );
    });
    gate.request();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1).unwrap();
        assert!(matches!(
            scaling.free_surface3_from_tetrahedra(domain(), [0.0; 3], snapshot, limits(), cx),
            Err(TetLiquidTransferError::Cancelled)
        ));
    });
}
