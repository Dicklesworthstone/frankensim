//! SI liquid impulses use actual partial-fill inertia and never integrate dt.
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_lbm::d3q19::freesurface3::ContactModel3;
use fs_lbm::freesurface::ContactModel;
use fs_lbm::impulse::ImpulseError;
use fs_lbm::material::free_surface::{
    MaterialFreeSurface2, MaterialFreeSurface3, MaterialFreeSurfaceDomain2,
    MaterialFreeSurfaceDomain3, MaterialFreeSurfaceImpulseControls2,
    MaterialFreeSurfaceImpulseControls3, MaterialFreeSurfaceImpulseError,
    plan_free_surface_from_material,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::surface::{
    LiquidCapillarity, ResolvedLiquidSurfaceState, resolve_liquid_surface_state,
};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::Dims;

fn liquid() -> ResolvedLiquidSurfaceState {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("density", Dims([-3, 1, 0, 0, 0, 0]), 1000.0),
        ("dynamic_viscosity", Dims([-1, 1, -1, 0, 0, 0]), 0.001),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 300.0)
                .with("surface_environment", 2.0, 2.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic liquid impulse fixture".into(),
                license: "test-only".into(), artifact: None,
            },
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic liquid".into(), phase: "liquid".into(),
        process: "numerical fixture".into(), revision: 0,
    }, claims, vec![]).unwrap();
    let point = QueryPoint::new().with("T", 300.0).unwrap()
        .with("surface_environment", 2.0).unwrap();
    resolve_liquid_surface_state(
        &card, &point, LiquidViscosityBasis::Dynamic, LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap()
}
fn flow3(source: &ResolvedLiquidSurfaceState, dt: f64) -> MaterialFreeSurface3<'_> {
    plan_free_surface_from_material(source, 0.001, dt, 0.1).unwrap()
        .free_surface3(MaterialFreeSurfaceDomain3 {
            cells: [4; 3], gravity_m_s2: [0.1, -0.2, 0.05],
            contact: ContactModel3::Neutral, solid_tiles: &[], max_cells: 64,
        }, |x, _, _| x < 2).unwrap()
}
fn flow2(source: &ResolvedLiquidSurfaceState, dt: f64) -> MaterialFreeSurface2<'_> {
    plan_free_surface_from_material(source, 0.001, dt, 0.1).unwrap()
        .free_surface2(MaterialFreeSurfaceDomain2 {
            cells: [6; 2], gravity_m_s2: [0.1, -0.2],
            contact: ContactModel::Neutral, max_cells: 36,
        }, |x, _| x < 3).unwrap()
}
fn controls3() -> MaterialFreeSurfaceImpulseControls3 {
    MaterialFreeSurfaceImpulseControls3 {
        max_cells: 64, maximum_speed_m_s: 0.08,
        maximum_momentum_defect_kg_m_s: 1e-16, maximum_work_defect_j: 1e-18,
    }
}
fn controls2() -> MaterialFreeSurfaceImpulseControls2 {
    MaterialFreeSurfaceImpulseControls2 {
        max_cells: 36, maximum_speed_m_s: 0.08,
        maximum_momentum_defect_kg_s: 1e-13, maximum_work_defect_j_m: 1e-15,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}");
}

#[test]
fn si_total_impulse_and_signed_work_use_liquid_mass_without_timestep_scaling() {
    let source = liquid();
    let gate = CancelGate::new_clock_free();
    for dt in [0.001, 0.0005] {
        let mut flow = flow3(&source, dt);
        let cells = [[0, 1, 1], [1, 1, 1]];
        let before = cells.map(|[x, y, z]| flow.sample_at_si(x, y, z).unwrap().unwrap());
        assert_eq!(before.map(|s| s.fill_fraction), [1.0, 0.5]);
        let mass = flow.ledger_mass_kg().unwrap().to_bits();
        let step = flow.lattice().grid().steps();
        let j = [10.0, -5.0, 2.0];
        let mut field = vec![[0.0; 3]; 64];
        field[20] = j;
        field[21] = j;
        let report = flow.apply_impulse_density(&field, controls3(), &gate).unwrap();
        let after = cells.map(|[x, y, z]| flow.sample_at_si(x, y, z).unwrap().unwrap());
        let mut observed = [0.0; 3];
        let mut work = 0.0;
        for (a, b) in before.into_iter().zip(after) {
            assert_eq!(a.cell, b.cell);
            assert_eq!(a.fill_fraction.to_bits(), b.fill_fraction.to_bits());
            assert_eq!(a.density_kg_m3.to_bits(), b.density_kg_m3.to_bits());
            for axis in 0..3 {
                let jump = b.velocity_m_s[axis] - a.velocity_m_s[axis];
                near(jump, j[axis] / (a.density_kg_m3 * a.fill_fraction), 5e-14);
                observed[axis] += a.density_kg_m3 * a.fill_fraction * 1e-9 * jump;
                work += j[axis] * 1e-9 * 0.5 * (a.velocity_m_s[axis] + b.velocity_m_s[axis]);
            }
        }
        for axis in 0..3 {
            near(report.applied_impulse_kg_m_s[axis], 2.0 * j[axis] * 1e-9, 1e-22);
            near(report.momentum_change_kg_m_s[axis], observed[axis], 1e-20);
        }
        near(report.work_j, work, 1e-20);
        near(report.kinetic_energy_change_j, work, 1e-20);
        assert_eq!(report.affected_cells, 2);
        assert_eq!(report.mass_change_kg, 0.0);
        assert_eq!(flow.ledger_mass_kg().unwrap().to_bits(), mass);
        assert_eq!(flow.lattice().grid().steps(), step);
        for impulse in &mut field { *impulse = impulse.map(|v| -v); }
        let braking = flow.apply_impulse_density(&field, controls3(), &gate).unwrap();
        assert!(braking.work_j < 0.0);
        near(report.work_j + braking.work_j, 0.0, 1e-20);
        assert!(std::ptr::eq(flow.scaling().material(), &source));
    }
}

#[test]
fn si_two_dimensional_impulse_and_work_use_original_units_per_depth() {
    let source = liquid();
    let gate = CancelGate::new_clock_free();
    for dt in [0.001, 0.0005] {
        let mut flow = flow2(&source, dt);
        let cells = [[1, 2], [2, 2]];
        let before = cells.map(|[x, y]| flow.sample_at_si(x, y).unwrap().unwrap());
        assert_eq!(before.map(|s| s.fill_fraction), [1.0, 0.5]);
        let mass = flow.ledger_mass_kg_per_m().unwrap().to_bits();
        let native_mass = flow.lattice().mass.clone();
        let flags = flow.lattice().grid.flags.clone();
        let force = (flow.lattice().grid.g, flow.lattice().grid.fext.clone());
        let j = [10.0, -5.0];
        let mut field = vec![[0.0; 2]; 36];
        field[13] = j;
        field[14] = j;
        let report = flow.apply_impulse_density(&field, controls2(), &gate).unwrap();
        let after = cells.map(|[x, y]| flow.sample_at_si(x, y).unwrap().unwrap());
        let mut observed = [0.0; 2];
        let mut work = 0.0;
        for (a, b) in before.into_iter().zip(after) {
            assert_eq!(a.fill_fraction.to_bits(), b.fill_fraction.to_bits());
            assert_eq!(a.density_kg_m3.to_bits(), b.density_kg_m3.to_bits());
            for axis in 0..2 {
                let jump = b.velocity_m_s[axis] - a.velocity_m_s[axis];
                near(jump, j[axis] / (a.density_kg_m3 * a.fill_fraction), 5e-14);
                observed[axis] += a.density_kg_m3 * a.fill_fraction * 1e-6 * jump;
                work += j[axis] * 1e-6 * 0.5 * (a.velocity_m_s[axis] + b.velocity_m_s[axis]);
            }
        }
        for axis in 0..2 {
            near(report.applied_impulse_kg_s[axis], 2.0 * j[axis] * 1e-6, 1e-19);
            near(report.momentum_change_kg_s[axis], observed[axis], 1e-17);
        }
        near(report.work_j_m, work, 1e-17);
        near(report.kinetic_energy_change_j_m, work, 1e-17);
        assert_eq!(report.affected_cells, 2);
        assert_eq!(report.mass_change_kg_m, 0.0);
        assert_eq!(flow.ledger_mass_kg_per_m().unwrap().to_bits(), mass);
        assert_eq!(flow.lattice().mass, native_mass);
        assert_eq!(flow.lattice().grid.flags, flags);
        assert_eq!(flow.lattice().grid.g, force.0);
        assert_eq!(flow.lattice().grid.fext, force.1);
        for impulse in &mut field { *impulse = impulse.map(|v| -v); }
        let braking = flow.apply_impulse_density(&field, controls2(), &gate).unwrap();
        assert!(braking.work_j_m < 0.0);
        near(report.work_j_m + braking.work_j_m, 0.0, 1e-17);
    }
}

#[test]
fn physical_conversion_scenario_bounds_and_cancellation_publish_nothing() {
    let source = liquid();
    let gate = CancelGate::new_clock_free();
    let mut three = flow3(&source, 0.001);
    let mut two = flow2(&source, 0.001);
    let original3 = three.lattice().grid().state_bits();
    let original2 = two.lattice().grid.f.clone();
    let mass3 = three.ledger_mass_kg().unwrap().to_bits();
    let mass2 = two.ledger_mass_kg_per_m().unwrap().to_bits();
    let mut field3 = vec![[0.0; 3]; 64];
    let mut field2 = vec![[0.0; 2]; 36];
    let mut wide3 = controls3();
    let mut wide2 = controls2();
    wide3.maximum_speed_m_s = 0.11;
    wide2.maximum_speed_m_s = 0.11;
    assert!(matches!(three.apply_impulse_density(&field3, wide3, &gate),
        Err(MaterialFreeSurfaceImpulseError::Units(_))));
    assert!(matches!(two.apply_impulse_density(&field2, wide2, &gate),
        Err(MaterialFreeSurfaceImpulseError::Units(_))));
    for unrepresentable in [f64::NAN, f64::from_bits(1)] {
        field3[20] = [10.0, 0.0, 0.0];
        field3[21][0] = unrepresentable;
        field2[13] = [10.0, 0.0];
        field2[14][0] = unrepresentable;
        assert!(matches!(three.apply_impulse_density(&field3, controls3(), &gate),
            Err(MaterialFreeSurfaceImpulseError::Units(_))));
        assert!(matches!(two.apply_impulse_density(&field2, controls2(), &gate),
            Err(MaterialFreeSurfaceImpulseError::Units(_))));
    }
    gate.request();
    assert_eq!(three.apply_impulse_density(&field3, controls3(), &gate),
        Err(MaterialFreeSurfaceImpulseError::Native(ImpulseError::Cancelled)));
    assert_eq!(two.apply_impulse_density(&field2, controls2(), &gate),
        Err(MaterialFreeSurfaceImpulseError::Native(ImpulseError::Cancelled)));
    assert_eq!(three.lattice().grid().state_bits(), original3);
    assert_eq!(two.lattice().grid.f, original2);
    assert_eq!(three.ledger_mass_kg().unwrap().to_bits(), mass3);
    assert_eq!(two.ledger_mass_kg_per_m().unwrap().to_bits(), mass2);
}
