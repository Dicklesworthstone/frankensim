//! G1/G3 physical pressure, gauge and density substitution in the dense owner.
use super::*;
use fs_exec::CancelGate;
use fs_lbm::core2::Cell;
use fs_lbm::freesurface::{ContactModel, FreeSurfacePressureError2};
use fs_lbm::material::MaterialLbmError;
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

fn pressure_domain(gravity_m_s2: [f64; 2]) -> MaterialFreeSurfaceDomain2 {
    MaterialFreeSurfaceDomain2 {
        cells: [12; 2],
        gravity_m_s2,
        contact: ContactModel::Neutral,
        max_cells: 144,
    }
}

#[test]
fn g1_g3_initial_dense_pressure_matches_loaded_rest_and_material_units() {
    let gate = CancelGate::new_clock_free();
    let mut baseline = None;
    for density in [1000.0, 2000.0] {
        let source = liquid(density, density * 0.0001, None);
        let load_pa = density * 0.002;
        let initialize = |reference| {
            plan(&source)
                .with_pressure_reference(reference)
                .unwrap()
                .free_surface2_with_pressure(
                    pressure_domain([0.0; 2]),
                    |_, y| y < 7,
                    0.05,
                    |_| reference + load_pa,
                    &gate,
                )
                .unwrap()
        };
        let mut loaded = initialize(100_000.0);
        let shifted = initialize(104_096.0);
        assert_eq!(loaded.lattice().grid.f, shifted.lattice().grid.f);
        assert_eq!(loaded.lattice().mass, shifted.lattice().mass);
        assert_eq!(loaded.lattice().grid.flags, shifted.lattice().grid.flags);
        if let Some(previous) = &baseline {
            assert_eq!(&loaded.lattice().grid.f, previous);
        }
        baseline = Some(loaded.lattice().grid.f.clone());

        let rho = 1.0 + 0.002 / CS2;
        // Fifty full cells and ten half-filled interface cells, per metre depth.
        let mass = 55.0 * rho * density * 1e-6;
        near(loaded.ledger_mass_kg_per_m().unwrap(), mass, 2e-13);
        for _ in 0..8 {
            let report = loaded
                .step_with_gas_pressure_and_momentum(&gate, 0.05, |_| 100_000.0 + load_pa)
                .unwrap();
            assert_eq!(report.native.routes, 0);
            near(loaded.ledger_mass_kg_per_m().unwrap(), mass, 2e-13);
            for y in 1..7 {
                for x in 1..11 {
                    let sample = loaded.sample_at_si(x, y).unwrap().unwrap();
                    near(sample.density_kg_m3, rho * density, 3e-13);
                    near(sample.pressure_perturbation_pa, load_pa, 1e-10);
                    assert!(sample.velocity_m_s.iter().all(|u| u.abs() < 3e-13));
                }
            }
        }
    }
}

#[test]
fn g1_dense_initial_pressure_balances_held_gravity_in_every_axis_and_sign() {
    let gate = CancelGate::new_clock_free();
    let source = liquid(1000.0, 0.1, None);
    let reference = 100_000.0;
    for axis in 0..2 {
        for sign in [-1.0_f64, 1.0] {
            let g = sign * 0.0002;
            let mut gravity = [0.0; 2];
            gravity[axis] = g / 0.001; // g_lattice = g_SI*dt^2/dx.
            let surface_index = if sign < 0.0 { 6 } else { 5 };
            let fluid = |x, y| {
                let coordinate = [x, y][axis];
                if sign < 0.0 { coordinate < 7 } else { coordinate >= 5 }
            };
            // Independent stationary solution of the dense rho*g owner:
            // rho[j+1]/rho[j]=(1+a)/(1-a), a=g/(2*cs^2).
            // The missing-link pair at the free boundary requires
            // rho_gas=rho_interface*(1-|a|). This numerical profile is not an EOS.
            let a = g / (2.0 * CS2);
            let ratio = (1.0 + a) / (1.0 - a);
            let rho_at = |coordinate: i32| {
                (1.0 - a.abs()).recip() * ratio.powi(coordinate - surface_index)
            };
            let mut flow = plan(&source)
                .with_pressure_reference(reference)
                .unwrap()
                .free_surface2_with_pressure(
                    pressure_domain(gravity),
                    fluid,
                    0.05,
                    |position| {
                        let coordinate = (position[axis] / 0.001 - 0.5).round() as i32;
                        reference + 1000.0 * CS2 * (rho_at(coordinate) - 1.0)
                    },
                    &gate,
                )
                .unwrap();
            let mass = flow.ledger_mass_kg_per_m().unwrap();
            for _ in 0..16 {
                let report = flow
                    .step_with_gas_pressure_and_momentum(&gate, 0.05, |_| reference)
                    .unwrap();
                assert_eq!(report.native.routes, 0);
                near(flow.ledger_mass_kg_per_m().unwrap(), mass, 2e-12);
                for y in 1..11 {
                    for x in 1..11 {
                        let i = flow.lattice().grid.idx(x, y);
                        if !matches!(flow.lattice().grid.flags[i], Cell::Fluid | Cell::Interface) {
                            continue;
                        }
                        let sample = flow.sample_at_si(x, y).unwrap().unwrap();
                        let expected = 1000.0 * rho_at([x, y][axis] as i32);
                        near(sample.density_kg_m3, expected, 3e-12);
                        assert!(sample.velocity_m_s.iter().all(|u| u.abs() < 5e-12));
                    }
                }
            }
        }
    }
}

#[test]
fn g3_zero_initial_pressure_replays_original_forced_capillary_construction() {
    let gate = CancelGate::new_clock_free();
    let source = liquid(1000.0, 0.1, Some(0.05));
    let scaling = plan(&source).with_pressure_reference(100_000.0).unwrap();
    let domain = pressure_domain([0.15, -0.25]);
    let ordinary = scaling.free_surface2(domain, |_, y| y < 7).unwrap();
    let initialized = scaling
        .free_surface2_with_pressure(domain, |_, y| y < 7, 0.05, |_| 100_000.0, &gate)
        .unwrap();
    assert_eq!(ordinary.lattice().grid.f, initialized.lattice().grid.f);
    assert_eq!(ordinary.lattice().grid.flags, initialized.lattice().grid.flags);
    assert_eq!(ordinary.lattice().mass, initialized.lattice().mass);
    assert_eq!(ordinary.lattice().grid.g, initialized.lattice().grid.g);
    assert_eq!(ordinary.lattice().grid.fext, initialized.lattice().grid.fext);
    assert_eq!(ordinary.lattice().sigma, initialized.lattice().sigma);
    assert_eq!(ordinary.lattice().carry, initialized.lattice().carry);
    assert_eq!(ordinary.ledger_mass_kg_per_m().unwrap(), initialized.ledger_mass_kg_per_m().unwrap());
}

#[test]
fn g0_g4_dense_pressure_initialization_refuses_missing_invalid_and_cancelled_state() {
    let source = liquid(1000.0, 0.1, Some(0.05));
    let gate = CancelGate::new_clock_free();
    let domain = pressure_domain([0.0; 2]);
    assert!(
        plan(&source)
            .free_surface2_with_pressure(domain, |_, y| y < 7, 0.05, |_| {
                panic!("missing datum must refuse before sampling");
            }, &gate)
            .is_err()
    );
    let scaling = plan(&source).with_pressure_reference(100_000.0).unwrap();
    for bound in [0.0, -1.0, 1.0, f64::NAN] {
        assert!(
            scaling
                .free_surface2_with_pressure(domain, |_, y| y < 7, bound, |_| {
                    panic!("invalid bound must refuse before sampling");
                }, &gate)
                .is_err()
        );
    }
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, 100_100.0] {
        let visits = std::cell::Cell::new(0);
        let result = scaling.free_surface2_with_pressure(
            domain,
            |_, y| y < 7,
            0.05,
            |_| {
                visits.set(visits.get() + 1);
                if visits.get() == 60 { bad } else { 100_001.0 }
            },
            &gate,
        );
        assert!(matches!(
            result,
            Err(MaterialLbmError::FreeSurfacePressure2(
                FreeSurfacePressureError2::Invalid { coordinates: Some([10, 6]), .. }
            ))
        ));
        assert_eq!(visits.get(), 60);
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        scaling.free_surface2_with_pressure(domain, |_, y| y < 7, 0.05, |_| {
            panic!("cancelled construction must not sample pressure");
        }, &cancelled),
        Err(MaterialLbmError::FreeSurfacePressure2(FreeSurfacePressureError2::Cancelled))
    ));
    let during = CancelGate::new_clock_free();
    let visits = std::cell::Cell::new(0);
    assert!(matches!(
        scaling.free_surface2_with_pressure(
            domain,
            |_, y| y < 7,
            0.05,
            |_| {
                visits.set(visits.get() + 1);
                if visits.get() == 25 { during.request(); }
                100_001.0
            },
            &during,
        ),
        Err(MaterialLbmError::FreeSurfacePressure2(FreeSurfacePressureError2::Cancelled))
    ));
    assert!((25..=30).contains(&visits.get()), "cancellation must stop at the next row");
}
