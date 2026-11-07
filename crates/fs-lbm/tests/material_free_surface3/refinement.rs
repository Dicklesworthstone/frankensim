//! Original material/native-owner continuation across physical resolution change.
//! Synthetic coefficients establish numerical coupling, not calibrated accuracy.
use super::*;
use fs_lbm::material::free_surface::{MaterialDyadicRefinementControls3, MaterialFreeSurface3};

fn controls() -> MaterialDyadicRefinementControls3 {
    MaterialDyadicRefinementControls3 {
        max_cells: 512,
        max_active_tiles: 8,
        maximum_speed_m_s: 0.08,
        maximum_mass_defect_kg: 1e-15,
        maximum_first_moment_defect_kg_m: 1e-17,
        maximum_momentum_defect_kg_m_s: 1e-16,
        maximum_angular_momentum_defect_kg_m2_s: 1e-18,
        maximum_kinetic_energy_defect_j: 1e-18,
        maximum_viscous_stress_defect_pa: 1e-8,
    }
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

// Independent SI observation: actual partial fill and intrinsic density,
// physical child centers, and physical velocities from the original owner.
fn observed(flow: &MaterialFreeSurface3<'_>) -> [f64; 11] {
    let [nx, ny, nz] = flow.lattice().domain_cells();
    let dx = flow.scaling().cell_size_m();
    let mut result = [0.0; 11];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let Some(sample) = flow.sample_at_si(x, y, z).unwrap() else {
                    continue;
                };
                let mass = sample.fill_fraction * sample.density_kg_m3 * dx.powi(3);
                let position = [x, y, z].map(|i| (i as f64 + 0.5) * dx);
                let p = sample.velocity_m_s.map(|u| mass * u);
                let angular = [
                    position[1] * p[2] - position[2] * p[1],
                    position[2] * p[0] - position[0] * p[2],
                    position[0] * p[1] - position[1] * p[0],
                ];
                result[0] += mass;
                for a in 0..3 {
                    result[1 + a] += mass * position[a];
                    result[4 + a] += p[a];
                    result[7 + a] += angular[a];
                }
                result[10] += 0.5 * mass * sample.velocity_m_s.iter().map(|u| u * u).sum::<f64>();
            }
        }
    }
    result
}

#[test]
fn evolved_liquid_refines_in_original_si_units_and_continues_native_evolution() {
    let source = liquid(0, 300.0, Some(0.05));
    let angle = wetting(&source, FRAC_PI_2);
    let mut coarse = plan(&source)
        .with_planar_wetting(PlanarWall3::ZMin, &angle)
        .unwrap()
        .free_surface3(
            MaterialFreeSurfaceDomain3 {
                cells: [4; 3],
                max_cells: 64,
                gravity_m_s2: [0.1, -0.2, 0.05],
                ..domain()
            },
            |x, _, _| x < 2,
        )
        .unwrap();
    let gate = CancelGate::new_clock_free();
    for _ in 0..3 {
        coarse.step_with_gate(&gate).unwrap();
    }
    let source_bits = coarse.lattice().grid().state_bits();
    let source_mass = coarse.ledger_mass_kg().unwrap().to_bits();
    let before = observed(&coarse);
    assert!(before[4..7].iter().any(|p| p.abs() > 1e-12));
    let (mut fine, report) = coarse.refine_dyadic(controls(), &gate).unwrap();
    let after = observed(&fine);
    assert_eq!(fine.lattice().domain_cells(), [8; 3]);
    assert!(std::ptr::eq(fine.scaling().material(), &source));
    assert_eq!(
        fine.scaling().planar_wetting(),
        coarse.scaling().planar_wetting()
    );
    assert_eq!(
        fine.scaling().cell_size_m(),
        0.5 * coarse.scaling().cell_size_m()
    );
    assert_eq!(
        fine.scaling().timestep_s(),
        0.5 * coarse.scaling().timestep_s()
    );
    assert_eq!(
        fine.scaling().surface_tension_lattice(),
        2.0 * coarse.scaling().surface_tension_lattice()
    );
    let viscosity = |flow: &MaterialFreeSurface3<'_>| {
        flow.scaling().flow().viscosity * flow.scaling().cell_size_m().powi(2)
            / flow.scaling().timestep_s()
    };
    near(viscosity(&fine), viscosity(&coarse), 1e-18);
    for i in 0..11 {
        near(after[i], before[i], if i == 0 { 1e-17 } else { 1e-19 });
    }
    near(report.before.mass_kg, before[0], 1e-17);
    near(report.after.mass_kg, after[0], 1e-17);
    for a in 0..3 {
        near(report.after.first_mass_moment_kg_m[a], after[1 + a], 1e-19);
        near(report.after.momentum_kg_m_s[a], after[4 + a], 1e-20);
        near(
            report.after.angular_momentum_kg_m2_s[a],
            after[7 + a],
            1e-22,
        );
    }
    for z in 0..4 {
        for y in 0..4 {
            for x in 0..4 {
                let Some(parent) = coarse.sample_at_si(x, y, z).unwrap() else {
                    continue;
                };
                for child in 0..8 {
                    let child = fine
                        .sample_at_si(
                            2 * x + (child & 1),
                            2 * y + ((child >> 1) & 1),
                            2 * z + (child >> 2),
                        )
                        .unwrap()
                        .unwrap();
                    assert_eq!(child.cell, parent.cell);
                    assert_eq!(
                        child.fill_fraction.to_bits(),
                        parent.fill_fraction.to_bits()
                    );
                    assert_eq!(
                        child.density_kg_m3.to_bits(),
                        parent.density_kg_m3.to_bits()
                    );
                    assert_eq!(
                        child.pressure_perturbation_pa.to_bits(),
                        parent.pressure_perturbation_pa.to_bits()
                    );
                    for a in 0..3 {
                        near(child.velocity_m_s[a], parent.velocity_m_s[a], 1e-13);
                    }
                }
            }
        }
    }
    assert_eq!(coarse.lattice().grid().state_bits(), source_bits);
    assert_eq!(coarse.ledger_mass_kg().unwrap().to_bits(), source_mass);
    // A returned candidate must remain an executable original native owner.
    // Two fine steps span one old dt, without asserting trajectory identity.
    let fine_bits = fine.lattice().grid().state_bits();
    fine.step_with_gate(&gate).unwrap();
    fine.step_with_gate(&gate).unwrap();
    assert_ne!(fine.lattice().grid().state_bits(), fine_bits);
    near(fine.ledger_mass_kg().unwrap(), before[0], 2e-15);
    assert_eq!(coarse.lattice().grid().state_bits(), source_bits);
}

#[test]
fn physical_limits_and_cancel_refuse_without_changing_the_original_owner() {
    let source = liquid(0, 300.0, None);
    let coarse = plan(&source)
        .free_surface3(
            MaterialFreeSurfaceDomain3 {
                cells: [4; 3],
                max_cells: 64,
                ..domain()
            },
            |x, _, _| x < 2,
        )
        .unwrap();
    let original = coarse.lattice().grid().state_bits();
    let gate = CancelGate::new_clock_free();
    let mut limited = controls();
    limited.max_cells = 511;
    assert!(coarse.refine_dyadic(limited, &gate).is_err());
    limited = controls();
    limited.maximum_speed_m_s = 0.11;
    assert!(coarse.refine_dyadic(limited, &gate).is_err());
    limited = controls();
    limited.maximum_angular_momentum_defect_kg_m2_s = f64::NAN;
    assert!(coarse.refine_dyadic(limited, &gate).is_err());
    gate.request();
    assert!(matches!(
        coarse.refine_dyadic(controls(), &gate),
        Err(MaterialLbmError::FreeSurface(FreeSurfaceError3::Cancelled))
    ));
    assert_eq!(coarse.lattice().grid().state_bits(), original);
    assert_eq!(coarse.scaling().cell_size_m(), 1e-3);
    assert_eq!(coarse.scaling().timestep_s(), 2e-4);
    assert!(std::ptr::eq(coarse.scaling().material(), &source));
}
