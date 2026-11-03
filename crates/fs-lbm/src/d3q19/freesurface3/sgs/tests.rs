use super::*;
use crate::d3q19::freesurface3::{ContactModel3, ConversionRedistribution3};
use crate::d3q19::{E3, Q3, W3, equilibrium3};
use crate::sgs3::SgsUse;

fn controls(cs: f64) -> SgsConfig {
    SgsConfig {
        model: SgsModel::Smagorinsky { cs },
        filter_width_lu: 1.7,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

fn assert_same(a: &FreeSurface3, b: &FreeSurface3) {
    assert_eq!(a.grid.state_bits(), b.grid.state_bits());
    assert_eq!(a.mass_snapshot(), b.mass_snapshot());
    assert_eq!(a.conversions(), b.conversions());
    assert_eq!(a.geometry_epoch(), b.geometry_epoch());
    assert_eq!(
        a.worst_step_violation().to_bits(),
        b.worst_step_violation().to_bits()
    );
    assert_eq!(a.local_force, b.local_force);
    for (&key, tile) in &a.tiles {
        assert_eq!(tile.cells, b.tiles[&key].cells);
        assert_eq!(tile.mass, b.tiles[&key].mass);
    }
}

// Independent Hermite population reconstruction with physical Guo half-source.
// The stress contains an authored isotropic part to detect dilatational
// contamination of the deviatoric SGS fit.
fn populations(
    rho: f64,
    u: [f64; 3],
    force: [f64; 3],
    tau: f64,
    strain: [[f64; 3]; 3],
    dilation: f64,
) -> [f64; Q3] {
    let mut f = equilibrium3(rho, u);
    for (q, e) in E3.iter().enumerate() {
        let e = [f64::from(e.0), f64::from(e.1), f64::from(e.2)];
        let eu: f64 = (0..3).map(|a| e[a] * u[a]).sum();
        let mut contraction = 0.0;
        for a in 0..3 {
            for b in 0..3 {
                let h = e[a] * e[b] - if a == b { 1.0 / 3.0 } else { 0.0 };
                let s = strain[a][b] + if a == b { dilation } else { 0.0 };
                contraction += h * (-2.0 * rho * tau * s / 3.0);
            }
        }
        let source: f64 = (0..3)
            .map(|a| (3.0 * (e[a] - u[a]) + 9.0 * eu * e[a]) * force[a])
            .sum();
        f[q] += W3[q] * (4.5 * contraction - 0.5 * source);
    }
    f
}

#[test]
fn g1_force_corrected_deviatoric_stress_solves_current_viscosity_and_drives_collision() {
    let gate = CancelGate::new_clock_free();
    let molecular_tau = 0.65;
    let molecular_nu = (molecular_tau - 0.5) / 3.0;
    let strain: [[f64; 3]; 3] = [
        [0.004, 0.002, -0.001],
        [0.002, -0.003, 0.0015],
        [-0.001, 0.0015, -0.001],
    ];
    let rate_squared: f64 = 2.0 * strain.iter().flatten().map(|s| s * s).sum::<f64>();
    let eddy = (0.2_f64 * 1.7).powi(2) * rate_squared.sqrt();
    let expected_tau = molecular_tau + 3.0 * eddy;
    for dilation in [0.0, 0.007] {
        let mut lattice = FreeSurface3::new(
            4,
            4,
            4,
            molecular_tau,
            [0.0; 3],
            0.0,
            ContactModel3::Neutral,
            |_, _, _| true,
        )
        .unwrap();
        let forces: BTreeMap<_, _> = (0..TILE_CELLS)
            .map(|lane| ((0, lane), [0.0004 + lane as f64 * 0.00001, -0.0003, 0.0002]))
            .collect();
        lattice.rebind_force_field(forces.clone(), &gate).unwrap();
        let mut expected_resolved = 0.0;
        let mut expected_modeled = 0.0;
        for lane in 0..TILE_CELLS {
            let rho = 0.7 + 0.05 * (lane % 4) as f64;
            let u = [0.06, -0.03, 0.02];
            lattice.grid.set_populations(
                0,
                lane,
                populations(rho, u, forces[&(0, lane)], expected_tau, strain, dilation),
            );
            expected_resolved += rho * molecular_nu * rate_squared;
            expected_modeled += rho * eddy * rate_squared;
        }
        let mut workspace = FreeSurfaceSgsWorkspace3::default();
        let report = lattice
            .prepare_sgs(&gate, controls(0.2), &mut workspace)
            .unwrap();
        assert_eq!(report.wet_cells, 64);
        assert_eq!(report.liquid_volume_lu, 64.0);
        assert!((report.maximum_eddy_viscosity_lu - eddy).abs() < 2e-14);
        assert!((report.resolved_dissipation_lu - expected_resolved).abs() < 1e-14);
        assert!((report.sgs_dissipation_lu - expected_modeled).abs() < 1e-14);
        for lane in 0..TILE_CELLS {
            let p = [lane % 4, (lane / 4) % 4, lane / 16];
            assert!((workspace.relaxation_at(p).unwrap() - expected_tau).abs() < 2e-13);
        }
        // Independently collide and pull every population into an interior cell.
        let mut expected = [0.0; Q3];
        for (q, e) in E3.iter().enumerate() {
            let p = [2 - e.0, 2 - e.1, 2 - e.2];
            let (slot, lane) = lattice
                .grid
                .resolve_source(i64::from(p[0]), i64::from(p[1]), i64::from(p[2]))
                .unwrap();
            let rho = 0.7 + 0.05 * (lane % 4) as f64;
            let u = [0.06, -0.03, 0.02];
            let f = lattice.grid.populations(slot, lane);
            let eq = equilibrium3(rho, u);
            let e = [f64::from(e.0), f64::from(e.1), f64::from(e.2)];
            let eu: f64 = (0..3).map(|a| e[a] * u[a]).sum();
            let source: f64 = (0..3)
                .map(|a| (3.0 * (e[a] - u[a]) + 9.0 * eu * e[a]) * forces[&(0, lane)][a])
                .sum();
            expected[q] =
                f[q] - (f[q] - eq[q]) / expected_tau + (1.0 - 0.5 / expected_tau) * W3[q] * source;
        }
        let mass = lattice.ledger_mass();
        lattice
            .step_with_sgs(&gate, controls(0.2), &mut workspace)
            .unwrap();
        let (slot, lane) = lattice.grid.resolve_source(2, 2, 2).unwrap();
        let actual = lattice.grid.populations(slot, lane);
        for q in 0..Q3 {
            assert!((actual[q] - expected[q]).abs() < 3e-14, "population {q}");
        }
        assert!((lattice.ledger_mass() - mass).abs() < 2e-12);
    }
}

#[test]
fn g0_g5_zero_coefficient_preserves_native_free_surface_trajectory_bitwise() {
    let gate = CancelGate::new_clock_free();
    for model in [
        ConversionRedistribution3::Local,
        ConversionRedistribution3::ConservativeLocal,
        ConversionRedistribution3::LegacyGlobal,
    ] {
        let mut ordinary = FreeSurface3::new(
            8,
            8,
            8,
            0.8,
            [1e-6, 0.0, -2e-6],
            0.001,
            ContactModel3::Neutral,
            |x, _, z| x < 5 && z < 5,
        )
        .unwrap()
        .with_conversion_redistribution(model);
        let mut zero = ordinary.clone();
        let mut workspace = FreeSurfaceSgsWorkspace3::default();
        for _ in 0..6 {
            ordinary.step_with_gate(&gate).unwrap();
            let report = zero
                .step_with_sgs(&gate, controls(0.0), &mut workspace)
                .unwrap();
            assert_eq!(report.sgs_dissipation_lu, 0.0);
            assert_eq!(report.maximum_eddy_viscosity_lu, 0.0);
            assert_same(&ordinary, &zero);
        }
    }
}

#[test]
fn g0_g5_active_sgs_moves_liquid_preserves_mass_and_replays() {
    let gate = CancelGate::new_clock_free();
    let mut initial = FreeSurface3::new(
        8,
        8,
        8,
        0.6,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        |x, _, z| x < 5 && z < 5,
    )
    .unwrap();
    for (&key, tile) in &initial.tiles {
        let slot = initial.grid.slot_of(key).unwrap();
        for lane in 0..TILE_CELLS {
            if tile.cells[lane] != Cell3::Gas {
                let (_, y, _) = FreeSurface3::global_of(key, lane);
                initial.grid.set_populations(
                    slot,
                    lane,
                    equilibrium3(1.0, [0.003 * (y as f64 - 3.5), 0.0, 0.015]),
                );
            }
        }
    }
    let mut a = initial.clone();
    let mut b = initial.clone();
    let mut ordinary = initial.clone();
    let mass = initial.ledger_mass();
    let mut work_a = FreeSurfaceSgsWorkspace3::default();
    let mut work_b = FreeSurfaceSgsWorkspace3::default();
    let mut modeled = 0.0;
    for _ in 0..10 {
        let ra = a.step_with_sgs(&gate, controls(0.2), &mut work_a).unwrap();
        let rb = b.step_with_sgs(&gate, controls(0.2), &mut work_b).unwrap();
        ordinary.step_with_gate(&gate).unwrap();
        modeled += ra.sgs_dissipation_lu;
        assert_eq!(ra, rb);
        assert_same(&a, &b);
        assert!((a.ledger_mass() - mass).abs() < 2e-10);
        assert!(ra.liquid_volume_lu < ra.wet_cells as f64);
    }
    assert!(modeled > 1e-10);
    assert_ne!(a.grid.state_bits(), ordinary.grid.state_bits());
}

#[test]
fn g0_g4_refusal_cancellation_and_sparse_activation_rebuild_only_scratch() {
    let gate = CancelGate::new_clock_free();
    let mut lattice = FreeSurface3::new(
        12,
        4,
        4,
        0.8,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        |x, _, _| x >= 8,
    )
    .unwrap();
    assert_eq!(lattice.retire_gas_tiles().unwrap(), 1);
    let (slot, lane) = lattice.grid.resolve_source(8, 1, 1).unwrap();
    let key = lattice.grid.active_keys()[slot];
    lattice.tiles.get_mut(&key).unwrap().mass[lane] = 1.05;
    lattice.grid.set_populations(
        slot,
        lane,
        populations(
            1.0,
            [0.0; 3],
            [0.0; 3],
            0.8,
            [[0.003, 0.002, 0.0], [0.002, -0.003, 0.0], [0.0; 3]],
            0.0,
        ),
    );
    let before = lattice.clone();
    let mass = lattice.ledger_mass();
    let mut workspace = FreeSurfaceSgsWorkspace3::default();
    for invalid in [
        SgsConfig {
            maximum_eddy_to_molecular: 0.0,
            ..controls(0.2)
        },
        SgsConfig {
            usage: SgsUse::AcceptanceRender,
            ..controls(0.2)
        },
        SgsConfig {
            model: SgsModel::Wale { cw: 0.5 },
            ..controls(0.2)
        },
        controls(f64::NAN),
    ] {
        assert!(
            lattice
                .step_with_sgs(&gate, invalid, &mut workspace)
                .is_err()
        );
        assert_same(&lattice, &before);
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        lattice.step_with_sgs(&cancelled, controls(0.2), &mut workspace),
        Err(FreeSurfaceSgsError3::Flow(FreeSurfaceError3::Cancelled))
    ));
    assert_same(&lattice, &before);
    lattice
        .step_with_sgs(&gate, controls(0.2), &mut workspace)
        .unwrap();
    assert!(lattice.conversions().tiles_activated > 0);
    assert!(lattice.conversions().gas_to_interface > 0);
    assert_ne!(lattice.grid.slot_of(key).unwrap(), slot);
    lattice
        .prepare_sgs(&gate, controls(0.2), &mut workspace)
        .unwrap();
    assert!(workspace.relaxation_at([7, 1, 1]).is_some());
    assert!(workspace.relaxation_at([0, 0, 0]).is_none());
    assert!((lattice.ledger_mass() - mass).abs() < 3e-12);
}
