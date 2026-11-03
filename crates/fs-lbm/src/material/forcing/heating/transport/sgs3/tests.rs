//! Native material/thermal coupling checks with manufactured source laws.
use super::*;
use crate::d3q19::equilibrium3;
use crate::material::{MaterialThermalSlab3, plan_scaling_from_material};
use crate::sgs3::dynamic::channel::{ChannelDynamicScalarConfig, ChannelDynamicSgsConfig};
use crate::sgs3::{SgsConfig, SgsError, SgsModel, SgsUse};
use std::cell::Cell;
use std::rc::Rc;

const PLANE: usize = 20;
const FLUID: usize = 100;
const CELLS: usize = FLUID + 2 * PLANE;

fn reference() -> FluidTransportState {
    FluidTransportState::try_new(2.0, 3.0, 0.02, 0.01).unwrap()
}
fn varying(t: f64) -> FluidTransportState {
    FluidTransportState::try_new(2.0, 3.0, 0.02 * (1.0 + 0.03 * (t - 300.0)),
        0.01 * (1.0 + 0.02 * (t - 300.0))).unwrap()
}
fn resolved(t: f64) -> FluidTransportState {
    let local = varying(t);
    FluidTransportState::try_new(2.0 * (1.0 - 0.002 * (t - 305.0)),
        3.0 * (1.0 + 0.001 * (t - 305.0)), local.thermal_conductivity_w_m_k(),
        local.dynamic_viscosity_pa_s()).unwrap()
}
fn local() -> ThermalSgsConfig {
    ThermalSgsConfig {
        closure: SgsConfig { model: SgsModel::Smagorinsky { cs: 0.15 },
            filter_width_lu: 1.0, maximum_eddy_to_molecular: 100.0, usage: SgsUse::Exploration },
        turbulent_prandtl: 0.9, maximum_mach: 0.2,
    }
}
fn channel() -> ThermalChannelSgsConfig {
    ThermalChannelSgsConfig {
        momentum: ChannelDynamicSgsConfig { maximum_coefficient_squared: 10.0,
            maximum_eddy_to_molecular: 100.0, usage: SgsUse::Exploration },
        heat: ChannelDynamicScalarConfig { maximum_coefficient: 10.0,
            maximum_diffusivity_to_molecular: 100.0, usage: SgsUse::Exploration },
        maximum_mach: 0.2,
    }
}
fn buoyancy() -> BoussinesqConfig<3> {
    BoussinesqConfig { gravity_m_s2: [0.0, 0.0, -9.81],
        max_relative_density_departure: 0.05, max_relative_specific_heat_departure: 0.05 }
}
fn limits() -> HeatSourceLimits {
    HeatSourceLimits { max_cells: CELLS, temperature_range_k: [290.0, 320.0],
        balance_tolerance_j: 1e-12 }
}
fn loads() -> (Vec<f64>, Vec<[f64; 3]>) {
    let mut power = vec![100.0; CELLS];
    power[..PLANE].fill(0.0);
    power[CELLS - PLANE..].fill(0.0);
    (power, vec![[0.1, 0.0, 0.03]; CELLS])
}
fn noise(i: usize, salt: u64) -> f64 {
    let mut h = (i as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ salt;
    h = (h ^ (h >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    h = (h ^ (h >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    h ^= h >> 31;
    (h >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
}
fn slab() -> MaterialSlab<ThermalLbm3> {
    let material = reference();
    let mut slab = plan_scaling_from_material(&material, 0.01, 1e-4, 0.1).unwrap()
        .thermal_slab3(MaterialThermalSlab3 { nx: 5, ny: 4, nz: 5,
            bottom_temperature_k: 310.0, top_temperature_k: 300.0,
            gravity_m_s2: 0.0, expansion_per_k: 0.001, max_cells: CELLS }).unwrap();
    for active in 0..FLUID {
        let cell = active + PLANE;
        let h = (active / PLANE) as f64 + 0.5;
        let envelope = 4.0 * h * (5.0 - h) / 25.0;
        let u = core::array::from_fn(|axis| 0.02 * envelope * noise(active, 31 * axis as u64 + 7));
        let accepted_force = [1e-5, -2e-5, 3e-5];
        slab.lattice.grid.f[cell] = equilibrium3(1.0, u);
        crate::thermal::forcing::reframe3(
            &mut slab.lattice.grid.f[cell], u, [0.0; 3], accepted_force,
        );
        slab.lattice.grid.fext[cell] = accepted_force;
        let t = 1.0 - h / 5.0 + 0.05 * envelope * noise(active, 101);
        let old: f64 = slab.lattice.gpop[cell].iter().sum();
        slab.lattice.gpop[cell][0] += t - old;
    }
    slab
}
fn bits(slab: &MaterialSlab<ThermalLbm3>) -> Vec<u64> {
    let s = &slab.lattice;
    s.grid.f.iter().flatten().chain(s.gpop.iter().flatten())
        .chain(s.grid.fext.iter().flatten()).chain(&s.grid.tau)
        .chain([&s.tau_g, &s.gbeta, &s.t_ref, &s.t_bottom, &s.t_top])
        .map(|x| x.to_bits()).collect()
}

#[derive(Debug)]
struct Ticket {
    value: FluidTransportState,
    query_k: f64,
    cell: usize,
    reads: Rc<Cell<usize>>,
}
impl FluidTransport for Ticket {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.reads.set(self.reads.get() + 1);
        Ok(self.value)
    }
}

#[test]
fn g0_material_sgs_constant_source_reproduces_both_native_closures() {
    let gate = CancelGate::new_clock_free();
    let (power, force) = loads();
    for dynamic in [false, true] {
        let mut actual = slab();
        let mut baseline = slab();
        let mut a = actual.lattice.thermal_sgs_workspace(CELLS).unwrap();
        let mut b = baseline.lattice.thermal_sgs_workspace(CELLS).unwrap();
        let mut ac = actual.lattice.thermal_channel_sgs_workspace(usize::MAX).unwrap();
        let mut bc = baseline.lattice.thermal_channel_sgs_workspace(usize::MAX).unwrap();
        for _ in 0..2 {
            let reads = Rc::new(Cell::new(0));
            let resolve = |cell, t, retained: &FluidTransportState| {
                assert_eq!(*retained, reference());
                Ok(Ticket { value: reference(), query_k: t, cell, reads: reads.clone() })
            };
            let (material, actual_report, expected_heat, expected_report) = if dynamic {
                let result = actual.step_with_material_transport_and_channel_sgs(
                    &power, &force, limits(), channel(), &mut ac, &gate, resolve).unwrap();
                let expected = baseline.step_with_power_density_and_channel_sgs(
                    &power, &force, limits(), channel(), &mut bc, &gate).unwrap();
                (result.material, result.sgs.transport, expected.heat, expected.sgs.transport)
            } else {
                let result = actual.step_with_material_transport_and_sgs(
                    &power, &force, limits(), local(), &mut a, &gate, resolve).unwrap();
                let expected = baseline.step_with_power_density_and_sgs(
                    &power, &force, limits(), local(), &mut b, &gate).unwrap();
                (result.material, result.sgs, expected.heat, expected.sgs)
            };
            assert_eq!(material.heat, expected_heat);
            assert_eq!(actual_report, expected_report);
            assert_eq!(bits(&actual), bits(&baseline));
            assert_eq!(reads.get(), 2 * FLUID);
            for sample in material.applied.iter().chain(&material.endpoint) {
                assert_eq!(sample.material.cell, sample.cell);
                assert_eq!(sample.material.query_k, sample.temperature_k);
                assert_eq!(sample.transport, reference());
            }
        }
    }
}

#[test]
fn g1_local_material_feedback_changes_transport_and_preserves_independent_channel_heat() {
    let gate = CancelGate::new_clock_free();
    let (power, force) = loads();
    for dynamic in [false, true] {
        let mut actual = slab();
        let mut frozen = slab();
        let original_rates = actual.lattice.grid.tau.clone();
        let original_thermal = actual.lattice.tau_g;
        let mut local_work = actual.lattice.thermal_sgs_workspace(CELLS).unwrap();
        let mut channel_work = actual.lattice.thermal_channel_sgs_workspace(usize::MAX).unwrap();
        let mut frozen_local = frozen.lattice.thermal_sgs_workspace(CELLS).unwrap();
        let mut frozen_channel = frozen.lattice.thermal_channel_sgs_workspace(usize::MAX).unwrap();
        for iteration in 0..2 {
            let original_t: Vec<_> = actual.lattice.gpop.iter()
                .map(|p| actual.temperature_si(p.iter().sum()).unwrap()).collect();
            let material = if dynamic {
                let result = actual.step_with_material_transport_and_channel_sgs(
                    &power, &force, limits(), channel(), &mut channel_work, &gate,
                    |_, t, _| Ok(varying(t))).unwrap();
                assert!(result.sgs.transport.maximum_eddy_diffusivity_lu > 0.0);
                if iteration == 0 {
                    let rate = channel_work.relaxation();
                    assert!(result.material.applied.iter().any(|s|
                        rate.momentum[s.cell] == s.momentum_relaxation
                        && rate.thermal[s.cell] > s.thermal_relaxation));
                }
                frozen.step_with_power_density_and_channel_sgs(
                    &power, &force, limits(), channel(), &mut frozen_channel, &gate).unwrap();
                result.material
            } else {
                let result = actual.step_with_material_transport_and_sgs(
                    &power, &force, limits(), local(), &mut local_work, &gate,
                    |_, t, _| Ok(varying(t))).unwrap();
                assert!(result.sgs.momentum.maximum_eddy_viscosity_lu > 0.0);
                frozen.step_with_power_density_and_sgs(
                    &power, &force, limits(), local(), &mut frozen_local, &gate).unwrap();
                result.material
            };
            for s in &material.applied {
                assert_eq!(s.temperature_k, original_t[s.cell]);
                assert_eq!(s.transport, varying(s.temperature_k));
            }
            for s in &material.endpoint {
                assert_eq!(s.temperature_k,
                    actual.temperature_si(actual.lattice.gpop[s.cell].iter().sum()).unwrap());
            }
            assert!(material.heat.defect_j.abs() <= limits().balance_tolerance_j);
            assert_eq!(actual.lattice.grid.tau, original_rates);
            assert_eq!(actual.lattice.tau_g, original_thermal);
        }
        assert_ne!(actual.lattice.grid.f, frozen.lattice.grid.f);
        assert_ne!(actual.lattice.gpop, frozen.lattice.gpop);
    }
    let mut actual = slab();
    let mut baseline = slab();
    let mut work = actual.lattice.thermal_sgs_workspace(CELLS).unwrap();
    let mut zero = local();
    zero.closure.model = SgsModel::Smagorinsky { cs: 0.0 };
    let result = actual.step_with_material_transport_and_sgs(
        &power, &force, limits(), zero, &mut work, &gate, |_, t, _| Ok(varying(t))).unwrap();
    let expected = baseline.step_with_material_transport(
        &power, &force, limits(), &gate, |_, t, _| Ok(varying(t))).unwrap();
    assert_eq!(result.material, expected);
    assert_eq!(bits(&actual), bits(&baseline));
}

#[test]
fn g1_buoyant_material_sgs_uses_source_force_once_and_retains_the_accepted_force_frame() {
    let gate = CancelGate::new_clock_free();
    let (power, external) = loads();
    for dynamic in [false, true] {
        let mut actual = slab();
        let mut reference = slab();
        actual.lattice.gbeta = 5e-5;
        let retained = actual.force_units;
        let mut projected_force = external.clone();
        for cell in PLANE..CELLS - PLANE {
            let t = actual.temperature_si(actual.lattice.gpop[cell].iter().sum()).unwrap();
            let departure = resolved(t).density_kg_m3() - 2.0;
            for axis in 0..3 {
                projected_force[cell][axis] += departure * buoyancy().gravity_m_s2[axis];
            }
        }
        let reads = Rc::new(Cell::new(0));
        let resolve = |cell, t, _: &FluidTransportState| {
            Ok(Ticket { value: resolved(t), query_k: t, cell, reads: reads.clone() })
        };
        let (material, report, expected_heat, expected_report) = if dynamic {
            let mut a = actual.lattice.thermal_channel_sgs_workspace(usize::MAX).unwrap();
            let mut b = reference.lattice.thermal_channel_sgs_workspace(usize::MAX).unwrap();
            let result = actual.step_with_boussinesq_material_transport_and_channel_sgs(
                &power, &external, limits(), buoyancy(), channel(), &mut a, &gate, resolve).unwrap();
            let expected = reference.step_with_material_transport_and_channel_sgs(
                &power, &projected_force, limits(), channel(), &mut b, &gate,
                |_, t, _| Ok(varying(t))).unwrap();
            (result.material, result.sgs.transport, expected.material.heat, expected.sgs.transport)
        } else {
            let mut a = actual.lattice.thermal_sgs_workspace(CELLS).unwrap();
            let mut b = reference.lattice.thermal_sgs_workspace(CELLS).unwrap();
            let result = actual.step_with_boussinesq_material_transport_and_sgs(
                &power, &external, limits(), buoyancy(), local(), &mut a, &gate, resolve).unwrap();
            let expected = reference.step_with_material_transport_and_sgs(
                &power, &projected_force, limits(), local(), &mut b, &gate,
                |_, t, _| Ok(varying(t))).unwrap();
            (result.material, result.sgs, expected.material.heat, expected.sgs)
        };
        assert_eq!(material.heat, expected_heat);
        assert_eq!(report, expected_report);
        assert_eq!(actual.lattice.grid.f, reference.lattice.grid.f);
        assert_eq!(actual.lattice.grid.fext, reference.lattice.grid.fext);
        assert_eq!(actual.lattice.gpop, reference.lattice.gpop);
        assert_eq!(actual.lattice.gbeta, 5e-5);
        assert!(actual.force_units == retained);
        assert_eq!(reads.get(), 2 * FLUID);
        for sample in material.applied.iter().chain(&material.endpoint) {
            assert_eq!(*sample.material.resolved_transport(), resolved(sample.temperature_k));
            assert_eq!(sample.material.material().query_k, sample.temperature_k);
            assert_eq!(sample.transport.density_kg_m3(), 2.0);
            assert_eq!(sample.transport.specific_heat_j_kg_k(), 3.0);
        }
    }
}

#[test]
fn g0_sgs_caps_use_fresh_molecular_source_properties_and_refusals_leave_no_step() {
    let gate = CancelGate::new_clock_free();
    let (power, force) = loads();
    let mut actual = slab();
    let before = bits(&actual);
    let mut work = actual.lattice.thermal_sgs_workspace(CELLS).unwrap();
    let fit = actual.lattice.prepare_thermal_sgs(local(), None, &mut work, &gate).unwrap();
    let mut config = local();
    config.closure.maximum_eddy_to_molecular = 2.0 * fit.momentum.maximum_eddy_to_molecular;
    let error = actual.step_with_material_transport_and_sgs(
        &power, &force, limits(), config, &mut work, &gate,
        |_, _, _| Ok(FluidTransportState::try_new(2.0, 3.0, 0.02, 0.001).unwrap())).unwrap_err();
    assert!(matches!(error, MaterialTransportError::Subgrid(
        ThermalSgsError::Closure(SgsError::ViscosityRatioExceeded { .. }))));
    assert_eq!(bits(&actual), before);
    let mut channel_work = actual.lattice.thermal_channel_sgs_workspace(usize::MAX).unwrap();
    let fit = actual.lattice.prepare_thermal_channel_sgs(
        channel(), None, &mut channel_work, &gate).unwrap();
    let mut config = channel();
    config.heat.maximum_diffusivity_to_molecular =
        2.0 * fit.transport.maximum_eddy_diffusivity_lu / ((actual.lattice.tau_g - 0.5) / 4.0);
    let error = actual.step_with_material_transport_and_channel_sgs(
        &power, &force, limits(), config, &mut channel_work, &gate,
        |_, _, _| Ok(FluidTransportState::try_new(2.0, 3.0, 0.002, 0.01).unwrap())).unwrap_err();
    assert!(matches!(error, MaterialTransportError::Subgrid(ThermalSgsError::DiffusivityRatio { .. })));
    assert_eq!(bits(&actual), before);
    let mut fresh = slab();
    let mut fresh_work = fresh.lattice.thermal_channel_sgs_workspace(usize::MAX).unwrap();
    let replay = actual.step_with_material_transport_and_channel_sgs(
        &power, &force, limits(), config, &mut channel_work, &gate, |_, _, _| Ok(reference())).unwrap();
    let expected = fresh.step_with_material_transport_and_channel_sgs(
        &power, &force, limits(), config, &mut fresh_work, &gate, |_, _, _| Ok(reference())).unwrap();
    assert_eq!(replay, expected);
    assert_eq!(bits(&actual), bits(&fresh));
}

#[test]
fn g0_g4_g5_material_endpoints_refuse_or_cancel_before_either_sgs_family_publishes() {
    let (power, force) = loads();
    for dynamic in [false, true] {
        for mode in 0..3 {
            let mut actual = slab();
            let before = bits(&actual);
            let gate = CancelGate::new_clock_free();
            let mut calls = 0;
            let mut resolve = |_: usize, t: f64, _: &FluidTransportState| {
                calls += 1;
                if calls > FLUID && mode == 0 {
                    return Err(MaterialError::State { what: "endpoint source refused".into() });
                }
                if calls > FLUID && mode == 1 {
                    return Ok(FluidTransportState::try_new(2.2, 3.0, 0.02, 0.01).unwrap());
                }
                if calls == 2 * FLUID && mode == 2 { gate.request(); }
                Ok(resolved(t))
            };
            let mut local_work = actual.lattice.thermal_sgs_workspace(CELLS).unwrap();
            let mut channel_work = actual.lattice.thermal_channel_sgs_workspace(usize::MAX).unwrap();
            let result = if dynamic {
                actual.step_with_boussinesq_material_transport_and_channel_sgs(
                    &power, &force, limits(), buoyancy(), channel(), &mut channel_work, &gate,
                    &mut resolve).map(|r| r.material.heat)
            } else {
                actual.step_with_boussinesq_material_transport_and_sgs(
                    &power, &force, limits(), buoyancy(), local(), &mut local_work, &gate,
                    &mut resolve).map(|r| r.material.heat)
            };
            if mode == 2 {
                assert!(matches!(result, Err(MaterialTransportError::Cancelled)));
                assert_eq!(calls, 2 * FLUID);
            } else {
                assert!(matches!(result, Err(MaterialTransportError::State { cell: PLANE, .. })));
                assert_eq!(calls, FLUID + 1);
            }
            assert_eq!(bits(&actual), before);
            let gate = CancelGate::new_clock_free();
            let mut fresh = slab();
            if dynamic {
                let mut fresh_work = fresh.lattice.thermal_channel_sgs_workspace(usize::MAX).unwrap();
                let replay = actual.step_with_boussinesq_material_transport_and_channel_sgs(
                    &power, &force, limits(), buoyancy(), channel(), &mut channel_work, &gate,
                    |_, t, _| Ok(resolved(t))).unwrap();
                let expected = fresh.step_with_boussinesq_material_transport_and_channel_sgs(
                    &power, &force, limits(), buoyancy(), channel(), &mut fresh_work, &gate,
                    |_, t, _| Ok(resolved(t))).unwrap();
                assert_eq!(replay, expected);
            } else {
                let mut fresh_work = fresh.lattice.thermal_sgs_workspace(CELLS).unwrap();
                let replay = actual.step_with_boussinesq_material_transport_and_sgs(
                    &power, &force, limits(), buoyancy(), local(), &mut local_work, &gate,
                    |_, t, _| Ok(resolved(t))).unwrap();
                let expected = fresh.step_with_boussinesq_material_transport_and_sgs(
                    &power, &force, limits(), buoyancy(), local(), &mut fresh_work, &gate,
                    |_, t, _| Ok(resolved(t))).unwrap();
                assert_eq!(replay, expected);
            }
            assert_eq!(bits(&actual), bits(&fresh));
        }
    }
}
