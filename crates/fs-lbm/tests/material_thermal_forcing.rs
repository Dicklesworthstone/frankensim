//! External multiphysics loads must survive the thermal buoyancy update.
use fs_lbm::ThermalLbm;
use fs_lbm::d3q19::ThermalLbm3;
use fs_lbm::material::{MaterialThermalSlab, MaterialThermalSlab3, plan_scaling_from_material};
use fs_material::fluid_transport::FluidTransportState;

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-10 * a.abs().max(b.abs()).max(1e-10), "{a} != {b}");
}

#[test]
fn external_forces_superpose_with_buoyancy_in_both_dimensions() {
    let mut two = ThermalLbm::slab(4, 6, 0.8, 0.7, 1e-5);
    let mut three = ThermalLbm3::slab(4, 3, 6, 0.8, 0.7, 1e-5);
    let i2 = two.grid.idx(1, 2);
    let i3 = three.grid.idx(1, 1, 2);
    let buoyancy2 = two.gbeta * (two.temperature(1, 2) - two.t_ref);
    let buoyancy3 = three.gbeta * (three.temperature(1, 1, 2) - three.t_ref);
    let loads2 = vec![[1e-6, -2e-6]; two.gpop.len()];
    let loads3 = vec![[1e-6, 3e-6, -2e-6]; three.gpop.len()];
    two.step_with_external_forces(&loads2).unwrap();
    three.step_with_external_forces(&loads3).unwrap();
    close(two.grid.fext[i2][0], 1e-6);
    close(two.grid.fext[i2][1], buoyancy2 - 2e-6);
    close(three.grid.fext[i3][0], 1e-6);
    close(three.grid.fext[i3][1], 3e-6);
    close(three.grid.fext[i3][2], buoyancy3 - 2e-6);
    assert!(two.grid.moments(i2).u[0] > 0.0);
    assert!(three.grid.moments3(i3).1[0] > 0.0);
}

#[test]
fn material_density_controls_real_force_response_without_callsite_conversion() {
    let mut velocities = Vec::new();
    for density in [1000.0, 2000.0] {
        let material = FluidTransportState::try_new(density, 4000.0, 0.6, 0.001).unwrap();
        let scale = plan_scaling_from_material(&material, 1e-4, 1e-3, 0.001).unwrap();
        let mut two = scale.thermal_slab(MaterialThermalSlab {
            nx: 4, ny: 6, bottom_temperature_k: 310.0, top_temperature_k: 300.0,
            gravity_m_s2: 0.0, expansion_per_k: 2e-4, max_cells: 32,
        }).unwrap();
        let mut three = scale.thermal_slab3(MaterialThermalSlab3 {
            nx: 4, ny: 3, nz: 6, bottom_temperature_k: 310.0, top_temperature_k: 300.0,
            gravity_m_s2: 0.0, expansion_per_k: 2e-4, max_cells: 96,
        }).unwrap();
        // Identical physical input; neither call site knows rho, nu or lattice units.
        let loads2 = vec![[100.0, 0.0]; two.lattice.gpop.len()];
        let loads3 = vec![[100.0, 0.0, 0.0]; three.lattice.gpop.len()];
        two.step_with_force_density(&loads2).unwrap();
        three.step_with_force_density(&loads3).unwrap();
        let i2 = two.lattice.grid.idx(1, 3);
        let i3 = three.lattice.grid.idx(1, 1, 3);
        close(two.lattice.grid.fext[i2][0], scale.force_density_lattice(100.0).unwrap());
        close(three.lattice.grid.fext[i3][0], scale.force_density_lattice(100.0).unwrap());
        let u2 = scale.velocity_si(two.lattice.grid.moments(i2).u[0]).unwrap();
        let u3 = scale.velocity_si(three.lattice.grid.moments3(i3).1[0]).unwrap();
        close(u2, u3);
        assert!(u2 > 0.0);
        velocities.push(u2);
    }
    close(velocities[0], 2.0 * velocities[1]);
}

#[test]
fn invalid_external_fields_leave_flow_heat_and_forces_unchanged() {
    let mut two = ThermalLbm::slab(4, 6, 0.8, 0.7, 1e-5);
    let mut three = ThermalLbm3::slab(4, 3, 6, 0.8, 0.7, 1e-5);
    let before2 = (two.grid.f.clone(), two.gpop.clone(), two.grid.fext.clone());
    let before3 = (three.grid.f.clone(), three.gpop.clone(), three.grid.fext.clone());
    assert!(two.step_with_external_forces(&[]).is_err());
    assert!(three.step_with_external_forces(&[]).is_err());
    let mut forces2 = vec![[0.0; 2]; two.gpop.len()];
    let mut forces3 = vec![[0.0; 3]; three.gpop.len()];
    forces2[5][0] = f64::NAN;
    forces3[15][2] = f64::INFINITY;
    assert!(two.step_with_external_forces(&forces2).is_err());
    assert!(three.step_with_external_forces(&forces3).is_err());
    forces2[5][0] = 0.0;
    forces3[15][2] = 0.0;
    // Finite inputs can still overflow when superposed with buoyancy.
    two.gbeta = f64::MAX; two.t_ref = -2.0;
    three.gbeta = f64::MAX; three.t_ref = -2.0;
    assert!(two.step_with_external_forces(&forces2).is_err());
    assert!(three.step_with_external_forces(&forces3).is_err());
    assert_eq!((two.grid.f, two.gpop, two.grid.fext), before2);
    assert_eq!((three.grid.f, three.gpop, three.grid.fext), before3);
}

#[test]
fn material_unit_refusals_do_not_partially_step_the_slab() {
    let mut slab = {
        let material = FluidTransportState::try_new(1000.0, 4000.0, 0.6, 0.001).unwrap();
        let scale = plan_scaling_from_material(&material, 1e-4, 1e-3, 0.001).unwrap();
        scale.thermal_slab3(MaterialThermalSlab3 {
            nx: 4, ny: 3, nz: 6, bottom_temperature_k: 310.0, top_temperature_k: 300.0,
            gravity_m_s2: 0.0, expansion_per_k: 2e-4, max_cells: 96,
        }).unwrap()
    }; // The slab owns its frozen units; a later plan cannot replace its density.
    let before = (slab.lattice.grid.f.clone(), slab.lattice.gpop.clone());
    assert!(slab.step_with_force_density(&[]).is_err());
    let mut loads = vec![[1.0, 0.0, 0.0]; slab.lattice.gpop.len()];
    loads[20][0] = f64::INFINITY;
    assert!(slab.step_with_force_density(&loads).is_err());
    assert_eq!((slab.lattice.grid.f.clone(), slab.lattice.gpop.clone()), before);
    loads[20][0] = 1.0;
    slab.step_with_force_density(&loads).unwrap();
    assert_ne!(slab.lattice.grid.f, before.0);
}

#[test]
fn explicit_zero_forcing_preserves_the_original_step_bits() {
    let mut a = ThermalLbm::slab(4, 6, 0.8, 0.7, 1e-5);
    let mut b = ThermalLbm::slab(4, 6, 0.8, 0.7, 1e-5);
    let mut c = ThermalLbm3::slab(4, 3, 6, 0.8, 0.7, 1e-5);
    let mut d = ThermalLbm3::slab(4, 3, 6, 0.8, 0.7, 1e-5);
    let zeros2 = vec![[0.0; 2]; b.gpop.len()];
    let zeros3 = vec![[0.0; 3]; d.gpop.len()];
    for _ in 0..5 {
        a.step(); b.step_with_external_forces(&zeros2).unwrap();
        c.step(); d.step_with_external_forces(&zeros3).unwrap();
    }
    assert!(a.grid.f.iter().flatten().map(|v| v.to_bits()).eq(b.grid.f.iter().flatten().map(|v| v.to_bits())));
    assert!(a.gpop.iter().flatten().map(|v| v.to_bits()).eq(b.gpop.iter().flatten().map(|v| v.to_bits())));
    assert!(c.grid.f.iter().flatten().map(|v| v.to_bits()).eq(d.grid.f.iter().flatten().map(|v| v.to_bits())));
    assert!(c.gpop.iter().flatten().map(|v| v.to_bits()).eq(d.gpop.iter().flatten().map(|v| v.to_bits())));
}

#[test]
fn external_loads_are_not_implicitly_accumulated_or_retained() {
    let mut two = ThermalLbm::slab(4, 6, 0.8, 0.7, 0.0);
    let mut three = ThermalLbm3::slab(4, 3, 6, 0.8, 0.7, 0.0);
    let forces2 = vec![[1e-6, 0.0]; two.gpop.len()];
    let forces3 = vec![[1e-6, 0.0, 0.0]; three.gpop.len()];
    for _ in 0..3 {
        two.step_with_external_forces(&forces2).unwrap();
        three.step_with_external_forces(&forces3).unwrap();
        assert_eq!(two.grid.fext[two.grid.idx(1, 3)][0], 1e-6);
        assert_eq!(three.grid.fext[three.grid.idx(1, 1, 3)][0], 1e-6);
    }
    two.step(); three.step();
    assert_eq!(two.grid.fext[two.grid.idx(1, 3)][0], 0.0);
    assert_eq!(three.grid.fext[three.grid.idx(1, 1, 3)][0], 0.0);
}

fn stress<const D: usize, const Q: usize>(
    populations: &[f64; Q],
    directions: [[f64; D]; Q],
    rho: f64,
    velocity: [f64; D],
    force: [f64; D],
) -> [[f64; D]; D] {
    core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            populations
                .iter()
                .zip(directions)
                .map(|(f, e)| f * e[a] * e[b])
                .sum::<f64>()
                + 0.5 * (velocity[a] * force[b] + force[a] * velocity[b])
                - rho * velocity[a] * velocity[b]
                - if a == b { rho / 3.0 } else { 0.0 }
        })
    })
}

#[test]
fn g1_buoyant_onset_perturbation_changes_only_the_requested_physical_velocity() {
    let mut two = ThermalLbm::slab(4, 6, 0.8, 0.7, 0.002);
    let mut three = ThermalLbm3::slab(4, 2, 6, 0.8, 0.7, 0.002);
    two.grid.g = [0.0003, -0.0002];
    let loads2 = vec![[0.0001, -0.0004]; two.gpop.len()];
    let loads3 = vec![[0.0001, 0.0003, -0.0004]; three.gpop.len()];
    for _ in 0..3 {
        two.step_with_external_forces(&loads2).unwrap();
        three.step_with_external_forces(&loads3).unwrap();
    }
    let bits2 = two.grid.f.clone();
    let bits3 = three.grid.f.clone();
    two.perturb(0.0);
    three.perturb(0.0);
    assert_eq!(two.grid.f, bits2);
    assert_eq!(three.grid.f, bits3);

    let i2 = two.grid.idx(1, 3);
    let i3 = three.grid.idx(1, 1, 3);
    // A nonzero shear defect must survive; replacing populations by equilibrium
    // would accidentally remove it even if the half-force velocity were fixed.
    two.grid.f[i2][1] += 0.001;
    two.grid.f[i2][3] += 0.001;
    two.grid.f[i2][0] -= 0.002;
    three.grid.f[i3][1] += 0.001;
    three.grid.f[i3][2] += 0.001;
    three.grid.f[i3][0] -= 0.002;
    let before2 = two.grid.moments(i2);
    let before3 = three.grid.moments3(i3);
    let force2 = core::array::from_fn(|a| two.grid.g[a].mul_add(before2.rho, two.grid.fext[i2][a]));
    let force3 = three.grid.fext[i3];
    let e2 = [
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [-1.0, 0.0],
        [0.0, -1.0],
        [1.0, 1.0],
        [-1.0, 1.0],
        [-1.0, -1.0],
        [1.0, -1.0],
    ];
    let e3 = fs_lbm::d3q19::E3.map(|(x, y, z)| [f64::from(x), f64::from(y), f64::from(z)]);
    let stress2 = stress(&two.grid.f[i2], e2, before2.rho, before2.u, force2);
    let stress3 = stress(&three.grid.f[i3], e3, before3.0, before3.1, force3);
    let amplitude = 0.001;
    let increment = amplitude
        * fs_math::det::sin(std::f64::consts::TAU / 4.0)
        * fs_math::det::sin(std::f64::consts::PI * 2.5 / 6.0);
    two.perturb(amplitude);
    three.perturb(amplitude);
    let after2 = two.grid.moments(i2);
    let after3 = three.grid.moments3(i3);
    assert!((after2.rho - before2.rho).abs() < 2e-14);
    assert!((after3.0 - before3.0).abs() < 2e-14);
    for a in 0..2 {
        let expected = before2.u[a] + if a == 1 { increment } else { 0.0 };
        assert!((after2.u[a] - expected).abs() < 2e-14);
    }
    for a in 0..3 {
        let expected = before3.1[a] + if a == 2 { increment } else { 0.0 };
        assert!((after3.1[a] - expected).abs() < 2e-14);
    }
    let after_stress2 = stress(&two.grid.f[i2], e2, after2.rho, after2.u, force2);
    let after_stress3 = stress(&three.grid.f[i3], e3, after3.0, after3.1, force3);
    for (a, b) in stress2.iter().flatten().zip(after_stress2.iter().flatten()) {
        assert!((a - b).abs() < 2e-14);
    }
    for (a, b) in stress3.iter().flatten().zip(after_stress3.iter().flatten()) {
        assert!((a - b).abs() < 2e-14);
    }
}
