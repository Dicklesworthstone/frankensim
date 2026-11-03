//! G1/G3 physical pressure, gauge and density substitution in the dense owner.
use super::*;
use fs_exec::CancelGate;
use fs_lbm::freesurface::ContactModel;
use fs_lbm::material::free_surface::{MaterialFreeSurface2, MaterialFreeSurfaceDomain2};

fn flow(source: &ResolvedLiquidSurfaceState, reference: Option<f64>) -> MaterialFreeSurface2<'_> {
    let plan = plan(source);
    let plan = match reference {
        Some(p) => plan.with_pressure_reference(p).unwrap(),
        None => plan,
    };
    plan.free_surface2(
        MaterialFreeSurfaceDomain2 {
            cells: [12; 2],
            gravity_m_s2: [0.0; 2],
            contact: ContactModel::Neutral,
            max_cells: 144,
        },
        |_, y| y < 7,
    )
    .unwrap()
}

#[test]
fn g1_g3_absolute_pressure_drives_dense_liquid_with_gauge_and_material_similarity() {
    let gate = CancelGate::new_clock_free();
    let mut baseline = None;
    for density in [1000.0, 2000.0] {
        let source = liquid(density, density * 0.0001, None);
        let mut plain = flow(&source, Some(100_000.0));
        let mut loaded = flow(&source, Some(100_000.0));
        let mut shifted = flow(&source, Some(104_096.0));
        let mass = loaded.ledger_mass_kg_per_m().unwrap();
        let pressure = density * 0.001;
        plain.step();
        loaded
            .step_with_gas_pressure(&gate, 0.05, |_| 100_000.0 + pressure)
            .unwrap();
        shifted
            .step_with_gas_link_pressure(&gate, 0.05, |link| {
                for d in 0..2 {
                    assert_eq!(
                        link.liquid_center_m[d],
                        (link.link.liquid_cell[d] as f64 + 0.5) * 0.001
                    );
                    assert_eq!(
                        link.gas_center_m[d],
                        (link.link.gas_cell[d] as f64 + 0.5) * 0.001
                    );
                    assert_eq!(
                        link.midpoint_m[d],
                        f64::midpoint(link.liquid_center_m[d], link.gas_center_m[d])
                    );
                }
                104_096.0 + pressure
            })
            .unwrap();
        assert_eq!(loaded.lattice().grid.f, shifted.lattice().grid.f);
        assert_eq!(loaded.lattice().mass, shifted.lattice().mass);
        assert_ne!(loaded.lattice().grid.f, plain.lattice().grid.f);
        assert!(loaded.sample_at_si(5, 6).unwrap().unwrap().velocity_m_s[1] < 0.0);
        near(loaded.ledger_mass_kg_per_m().unwrap(), mass, 2e-14);
        if let Some(previous) = &baseline {
            assert_eq!(&loaded.lattice().grid.f, previous);
        }
        baseline = Some(loaded.lattice().grid.f.clone());
    }
}

#[test]
fn g0_absolute_pressure_requires_reference_positive_load_and_atomic_admission() {
    let gate = CancelGate::new_clock_free();
    let source = liquid(1000.0, 0.1, Some(0.05));
    let mut missing = flow(&source, None);
    assert!(
        missing
            .step_with_gas_pressure(&gate, 0.05, |_| 100_000.0)
            .is_err()
    );
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, 101_000.0] {
        let mut state = flow(&source, Some(100_000.0));
        let before = state.lattice().clone();
        assert!(
            state
                .step_with_gas_link_pressure(&gate, 0.05, |link| {
                    if link.link.liquid_cell[0] == 10 {
                        bad
                    } else {
                        100_001.0
                    }
                })
                .is_err()
        );
        assert_eq!(state.lattice().grid.f, before.grid.f);
        assert_eq!(state.lattice().grid.flags, before.grid.flags);
        assert_eq!(state.lattice().mass, before.mass);
        assert_eq!(state.lattice().carry, before.carry);
    }
}
