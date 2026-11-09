//! Native dispersive shell refraction, coherent fields, and physical work.
use fs_bem::helmholtz::bc::{IncidentField, Monopole, PlaneWave};
use fs_bem::helmholtz::molecular::MolecularNearFieldOptions;
use fs_bem::helmholtz::transmission::molecular::MolecularFluidInterface;
use fs_bem::helmholtz::transmission::nested::{
    CoreShellFluidInterface, CoreShellInterface, CoreShellRegion,
    MolecularCoreShellFluidInterface,
};
use fs_bem::helmholtz::transmission::Options;
use fs_bem::helmholtz::HelmholtzError;
use fs_bem::{near_field, panel3d::SpherePanels};
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{RelaxationState, RelaxingMedium, input::read_relaxation};
use fs_math::c64::C64;

#[path = "support/molecular_core_shell_oracle.rs"]
mod oracle;
use oracle::MolecularShellOracle;

const OMEGA: f64 = 200.0;
fn state(temperature: f64, pressure: f64, relaxation_ms: Option<f64>) -> RelaxationState {
    // Manufactured caloric reservoir, not a measured air attenuation fit.
    let text = relaxation_ms.map_or_else(|| "declared 1 synthetic 0".into(), |ms| format!(
        "declared 1 synthetic 1 vibration 1 100 1 power {} 293.15 0 200 500", 101.325 * ms,
    ));
    let law = read_relaxation(&mut text.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, pressure).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    law.evaluate(&medium).unwrap()
}
fn plane() -> IncidentField<'static> {
    IncidentField::PlaneWave(PlaneWave { direction: [0.0, 0.0, 1.0], amplitude_pa: C64::ONE })
}
fn expik(k: C64, r: f64) -> C64 {
    C64::new((k.re * r).cos(), (k.re * r).sin()).scale((-k.im * r).exp())
}
fn close(actual: &[C64], expected: &[C64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    let scale = expected.iter().map(|z| z.abs()).fold(1e-30_f64, f64::max);
    for (&a, &b) in actual.iter().zip(expected) {
        assert!((a - b).abs() < tolerance * scale, "{a:?} != {b:?}");
    }
}

#[test]
fn molecular_concentric_shell_converges_to_independent_complex_waves_and_work() {
    let sources = [state(293.15, 101325.0, Some(0.1)), state(250.0, 101325.0, Some(5.0)), state(353.15, 101325.0, Some(2.0))];
    let k = sources.each_ref().map(|s| s.relaxing_wavenumber(OMEGA).unwrap());
    let rho = sources.each_ref().map(|s| s.gas().density);
    assert!(k.iter().all(|wave| wave.im > 0.0));
    let reference = MolecularShellOracle::new(OMEGA, [1.0, 0.45], rho, k, 16);
    let lower_order = MolecularShellOracle::new(OMEGA, [1.0, 0.45], rho, k, 14);
    assert!(reference.continuity_defect([1.0, 0.45]) < 2e-13);
    let exact_core = -reference.outward_work(2, 0.45);
    let exact_total = -reference.outward_work(1, 1.0);
    let exact_shell = exact_total - exact_core;
    assert!(exact_core > 0.0 && exact_shell > 0.0);
    let mut errors = Vec::new();
    let gate = CancelGate::new();
    for level in [1, 2] {
        let outer = SpherePanels::icosphere(1.0, level).unwrap();
        let core = SpherePanels::icosphere(0.45, level).unwrap();
        let owner = MolecularCoreShellFluidInterface::new(&outer, &core, sources.each_ref(), OMEGA, 0.01, Options::default(), &gate).unwrap();
        assert_eq!(owner.wavenumbers_per_m(), k);
        assert_eq!(owner.active_interface_count(), 2);
        for (actual, expected) in owner.material_sources().into_iter().zip(&sources) {
            assert!(std::ptr::eq(actual, expected));
        }
        let field = owner.scatter_field(plane(), 0.01, 100_000, &gate).unwrap();
        let mut mesh_errors = Vec::new();
        for (region, index, radius) in [(CoreShellRegion::Exterior, 0, 1.8), (CoreShellRegion::Shell, 1, 0.65), (CoreShellRegion::Core, 2, 0.15)] {
            let observed = field.pressure_velocity_at(region, &[[0.0, 0.0, radius]], 0.01, MolecularNearFieldOptions::default(), &gate).unwrap();
            let (p, v) = reference.axial(index, radius);
            let (pl, vl) = lower_order.axial(index, radius);
            assert!((p - pl).abs() < 2e-13 && (v - vl).abs() < 2e-15);
            mesh_errors.push((observed.scalar.pressure[0] - p).abs() / p.abs());
            mesh_errors.push((observed.particle_velocity_m_s[0][2] - v).abs() / v.abs());
            if index == 0 {
                let incident = expik(k[0], radius);
                let exact_scattered = p - incident;
                assert!(exact_scattered.abs() > 1e-4);
                mesh_errors.push((observed.scalar.pressure[0] - p).abs() / exact_scattered.abs());
            }
        }
        let actual_core = field.enclosed_absorbed_power_w(CoreShellInterface::Core).unwrap();
        let actual_shell = field.shell_absorbed_power_w().unwrap();
        if level == 2 {
            assert!(actual_core > 0.0 && actual_shell > 0.0);
        }
        let total = field.enclosed_absorbed_power_w(CoreShellInterface::Outer).unwrap();
        assert!((actual_core + actual_shell - total).abs() < 2e-16 * total.abs().max(1.0));
        mesh_errors.push((actual_core - exact_core).abs() / exact_core);
        mesh_errors.push((actual_shell - exact_shell).abs() / exact_shell);
        errors.push(mesh_errors);
    }
    for (quantity, (&coarse, &fine)) in errors[0].iter().zip(&errors[1]).enumerate() {
        assert!(fine < coarse, "quantity {quantity}: {coarse} -> {fine}");
        assert!(fine < 0.15, "quantity {quantity}: mesh error {fine}");
    }
}

#[test]
fn empty_reservoirs_recover_native_nested_fields_and_work() {
    let sources = [state(293.15, 101325.0, None), state(250.0, 101325.0, None), state(353.15, 101325.0, None)];
    let (outer, core) = (SpherePanels::icosphere(1.0, 0).unwrap(), SpherePanels::icosphere(0.45, 0).unwrap());
    let gate = CancelGate::new();
    let real = CoreShellFluidInterface::new(&outer, &core, [sources[0].gas(), sources[1].gas(), sources[2].gas()], OMEGA, 0.01, Options::default(), &gate).unwrap();
    let molecular = MolecularCoreShellFluidInterface::new(&outer, &core, sources.each_ref(), OMEGA, 0.01, Options::default(), &gate).unwrap();
    let a = real.scatter_field(plane(), 0.01, &gate).unwrap();
    let b = molecular.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    for boundary in [CoreShellInterface::Outer, CoreShellInterface::Core] {
        let (ta, tb) = (a.interface_trace(boundary).unwrap(), b.interface_trace(boundary).unwrap());
        close(tb.pressure_pa(), ta.pressure_pa(), 3e-12);
        close(tb.normal_velocity_m_s(), ta.normal_velocity_m_s(), 3e-11);
        assert!((tb.interface_work().surface_power - ta.interface_work().surface_power).abs() < 1e-13);
    }
    for (region, radius) in [(CoreShellRegion::Exterior, 1.8), (CoreShellRegion::Shell, 0.65), (CoreShellRegion::Core, 0.15)] {
        let points = [[0.0, 0.0, radius]];
        let va = a.pressure_velocity_at(region, &points, 0.01, near_field::Options::default(), &gate).unwrap();
        let vb = b.pressure_velocity_at(region, &points, 0.01, MolecularNearFieldOptions::default(), &gate).unwrap();
        close(&vb.scalar.pressure, &va.scalar.pressure, 2e-9);
        close(&vb.particle_velocity_m_s[0], &va.particle_velocity_m_s[0], 2e-9);
    }
}

#[test]
fn molecular_transparency_removes_only_matched_interfaces_and_keeps_original_fields() {
    let outside = state(293.15, 101325.0, Some(5.0));
    let inside = state(293.15, 101325.0, Some(0.5));
    let (outer, core) = (SpherePanels::icosphere(1.0, 0).unwrap(), SpherePanels::icosphere(0.45, 0).unwrap());
    let gate = CancelGate::new();
    // Equal equilibrium rho/c does NOT remove different reservoir responses.
    for transparent_outer in [true, false] {
        let sources = if transparent_outer { [&outside, &outside, &inside] } else { [&outside, &inside, &inside] };
        let (surface, active, removed) = if transparent_outer { (&core, CoreShellInterface::Core, CoreShellInterface::Outer) } else { (&outer, CoreShellInterface::Outer, CoreShellInterface::Core) };
        let nested = MolecularCoreShellFluidInterface::new(&outer, &core, sources, OMEGA, 0.01, Options::default(), &gate).unwrap();
        let single = MolecularFluidInterface::new(surface, &outside, &inside, OMEGA, Options::default(), &gate).unwrap();
        assert_eq!(nested.active_interface_count(), 1);
        let a = nested.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
        let b = single.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
        assert!(std::ptr::eq(a.radiating_surface(), surface));
        assert!(a.interface_trace(removed).is_none());
        assert!(a.shell_absorbed_power_w().is_none());
        close(a.scattered_pressure_pa(), b.scattered_pressure_pa(), 3e-11);
        close(a.scattered_normal_velocity_m_s(), b.scattered_normal_velocity_m_s(), 3e-11);
        close(a.interface_trace(active).unwrap().pressure_pa(), b.total_pressure_pa(), 3e-12);
        let points = [[0.0, 0.0, 0.15]];
        let va = a.pressure_velocity_at(CoreShellRegion::Core, &points, 0.01, MolecularNearFieldOptions::default(), &gate).unwrap();
        let vb = b.interior_pressure_velocity_adaptive_at(&points, 0.01, 1000, MolecularNearFieldOptions::default(), &gate).unwrap();
        close(&va.scalar.pressure, &vb.scalar.pressure, 3e-11);
        close(&va.particle_velocity_m_s[0], &vb.particle_velocity_m_s[0], 3e-11);
    }
    let owner = MolecularCoreShellFluidInterface::new(&outer, &core, [&outside; 3], OMEGA, 0.01, Options::default(), &gate).unwrap();
    assert_eq!(owner.active_interface_count(), 0);
    let field = owner.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    assert!(field.scattered_pressure_pa().iter().all(|&p| p == C64::ZERO));
    let k = owner.wavenumbers_per_m()[0];
    for (region, radius) in [(CoreShellRegion::Exterior, 1.8), (CoreShellRegion::Shell, 0.65), (CoreShellRegion::Core, 0.15)] {
        let observed = field.pressure_velocity_at(region, &[[0.0, 0.0, radius]], 0.01, MolecularNearFieldOptions::default(), &gate).unwrap();
        let p = expik(k, radius);
        close(&observed.scalar.pressure, &[p], 2e-15);
        close(&observed.particle_velocity_m_s[0], &[C64::ZERO, C64::ZERO, (k * p).scale(1.0 / (outside.gas().density * OMEGA))], 2e-15);
        assert_eq!(observed.scalar.quadrature_error_estimate_pa, vec![0.0]);
    }
}

#[test]
fn molecular_nested_refuses_invalid_regions_sources_pressure_and_work_budgets() {
    let outside = state(293.15, 101325.0, Some(5.0));
    let inside = state(353.15, 101325.0, Some(2.0));
    let wrong_pressure = state(293.15, 90000.0, Some(5.0));
    let (outer, core) = (SpherePanels::icosphere(1.0, 0).unwrap(), SpherePanels::icosphere(0.45, 0).unwrap());
    let gate = CancelGate::new();
    assert!(MolecularCoreShellFluidInterface::new(&outer, &core, [&outside, &wrong_pressure, &inside], OMEGA, 0.01, Options::default(), &gate).is_err());
    assert!(MolecularCoreShellFluidInterface::new(&outer, &core, [&outside, &inside, &outside], OMEGA, 0.01, Options { maximum_bytes: 1, ..Options::default() }, &gate).is_err());
    let owner = MolecularCoreShellFluidInterface::new(&outer, &core, [&outside, &inside, &outside], OMEGA, 0.01, Options::default(), &gate).unwrap();
    assert!(owner.scatter_field(plane(), 0.01, 39, &gate).is_err()); // 20 + 20 source points
    let emitters = [Monopole { position_m: [0.0; 3], volume_velocity_m3_s: C64::ONE }];
    let internal = IncidentField::Monopoles(&emitters);
    assert!(owner.scatter_field(internal, 0.01, 1000, &gate).is_err());
    let field = owner.scatter_field(plane(), 0.01, 1000, &gate).unwrap();
    assert!(field.pressure_velocity_at(CoreShellRegion::Shell, &[[0.0; 3]], 0.01, MolecularNearFieldOptions::default(), &gate).is_err());
    assert!(field.pressure_velocity_at(CoreShellRegion::Core, &[[0.0, 0.0, 1.8]], 0.01, MolecularNearFieldOptions::default(), &gate).is_err());
    let too_small = MolecularNearFieldOptions { maximum_row_bytes: 1, ..MolecularNearFieldOptions::default() };
    assert!(field.pressure_velocity_at(CoreShellRegion::Shell, &[[0.0, 0.0, 0.65]], 0.01, too_small, &gate).is_err());
    let one_surface = MolecularNearFieldOptions { quadrature: near_field::Options { maximum_kernel_evaluations: 80 * outer.areas().len(), ..near_field::Options::default() }, ..MolecularNearFieldOptions::default() };
    assert!(field.pressure_velocity_at(CoreShellRegion::Shell, &[[0.0, 0.0, 0.65]], 0.01, one_surface, &gate).is_err());
    gate.request();
    assert!(matches!(field.pressure_velocity_at(CoreShellRegion::Exterior, &[[0.0, 0.0, 1.8]], 0.01, MolecularNearFieldOptions::default(), &gate), Err(HelmholtzError::Cancelled)));
}
