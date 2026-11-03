//! G0/G1 identities on the actual source-coupled slab, not a turbulent DNS claim.

use super::*;
use crate::d3q19::equilibrium3;
use crate::sgs3::SgsUse;

fn config() -> ThermalChannelSgsConfig {
    ThermalChannelSgsConfig {
        momentum: ChannelDynamicSgsConfig {
            maximum_coefficient_squared: 10.0,
            maximum_eddy_to_molecular: 100.0,
            usage: SgsUse::Exploration,
        },
        heat: ChannelDynamicScalarConfig {
            maximum_coefficient: 10.0,
            maximum_diffusivity_to_molecular: 100.0,
            usage: SgsUse::Exploration,
        },
        maximum_mach: 0.2,
    }
}

fn limits(n: usize) -> ThermalSourceLimits {
    ThermalSourceLimits {
        max_cells: n,
        temperature_range: [-0.2, 1.2],
        balance_tolerance: 2e-12,
    }
}

fn bits(state: &ThermalLbm3) -> Vec<u64> {
    state
        .grid
        .f
        .iter()
        .flatten()
        .chain(state.gpop.iter().flatten())
        .chain(state.grid.fext.iter().flatten())
        .chain(&state.grid.tau)
        .chain([
            &state.tau_g,
            &state.t_bottom,
            &state.t_top,
            &state.gbeta,
            &state.t_ref,
        ])
        .map(|x| x.to_bits())
        .collect()
}

fn random(i: usize, salt: u64) -> f64 {
    let mut h = (i as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ salt;
    h ^= h >> 30;
    h = h.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    h ^= h >> 27;
    h = h.wrapping_mul(0x94d0_49bb_1331_11eb);
    h ^= h >> 31;
    (h >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
}

// An intentionally small broad-band field with positive and negative plane
// fits. It obeys the no-slip envelope but is not a turbulent inlet generator.
fn seed(state: &mut ThermalLbm3) {
    let (nx, ny, nz) = (state.grid.nx, state.grid.ny, state.grid.nz - 2);
    for z in 0..nz {
        let h = z as f64 + 0.5;
        let wall = 4.0 * h * (nz as f64 - h) / (nz * nz) as f64;
        for y in 0..ny {
            for x in 0..nx {
                let active = (z * ny + y) * nx + x;
                let cell = state.grid.idx(x, y, z + 1);
                let u = core::array::from_fn(|c| 0.02 * wall * random(active, 31 * c as u64 + 7));
                let force = [1e-5, -2e-5, 3e-5];
                state.grid.f[cell] = equilibrium3(1.0, u);
                crate::thermal::forcing::reframe3(&mut state.grid.f[cell], u, [0.0; 3], force);
                state.grid.fext[cell] = force;
                let temperature = 1.0 - h / nz as f64 + 0.05 * wall * random(active, 101);
                let old: f64 = state.gpop[cell].iter().sum();
                state.gpop[cell][0] += temperature - old;
            }
        }
    }
}

#[test]
fn g0_zero_channel_fit_retains_original_source_steps_and_local_rates_bitwise() {
    let gate = CancelGate::new_clock_free();
    for local in [false, true] {
        let mut reference = ThermalLbm3::slab(4, 3, 4, 0.8, 0.9, 0.0);
        let mut actual = reference.snapshot().0;
        let n = actual.gpop.len();
        let mut work = actual.thermal_channel_sgs_workspace(usize::MAX).unwrap();
        let pointer = work.relaxation().thermal.as_ptr();
        let momentum: Vec<_> = (0..n).map(|i| 0.7 + 0.01 * (i % 5) as f64).collect();
        let thermal: Vec<_> = (0..n).map(|i| 0.8 + 0.01 * (i % 7) as f64).collect();
        let molecular = local.then_some(RelaxationField {
            momentum: &momentum,
            thermal: &thermal,
        });
        let mut heat = vec![0.0001; n];
        heat[..12].fill(0.0);
        heat[n - 12..].fill(0.0);
        let force = vec![[0.0; 3]; n];
        let expected = match molecular {
            Some(field) => {
                reference.step_with_sources_and_transport(&heat, &force, field, limits(n))
            }
            None => reference.step_with_sources(&heat, &force, limits(n)),
        }
        .unwrap();
        let result = actual
            .step_with_sources_and_channel_sgs(
                &heat,
                &force,
                config(),
                molecular,
                &mut work,
                limits(n),
                &gate,
            )
            .unwrap();
        assert_eq!(result.heat, expected);
        assert_eq!(result.sgs.transport.momentum.maximum_eddy_viscosity_lu, 0.0);
        assert_eq!(result.sgs.transport.maximum_eddy_diffusivity_lu, 0.0);
        assert_eq!(bits(&actual), bits(&reference));
        assert_eq!(pointer, work.relaxation().thermal.as_ptr());
    }
}

#[test]
fn g1_channel_heat_survives_zero_momentum_and_reaches_original_transport() {
    let gate = CancelGate::new_clock_free();
    let mut state = ThermalLbm3::slab(5, 4, 5, 0.8, 0.9, 1e-5);
    seed(&mut state);
    let n = state.gpop.len();
    let original_rates = state.grid.tau.clone();
    let mut oracle = state.snapshot().0;
    let mut work = state.thermal_channel_sgs_workspace(usize::MAX).unwrap();
    let mut momentum = vec![0.8; n];
    let mut thermal = vec![0.9; n];
    for cell in 20..n - 20 {
        momentum[cell] += 0.001 * (cell % 5) as f64;
        thermal[cell] += 0.001 * (cell % 7) as f64;
    }
    let molecular = RelaxationField {
        momentum: &momentum,
        thermal: &thermal,
    };
    let prepared = state
        .prepare_thermal_channel_sgs(config(), Some(molecular), &mut work, &gate)
        .unwrap();
    assert!(prepared.transport.momentum.maximum_eddy_viscosity_lu > 0.0);
    assert!(prepared.transport.maximum_eddy_diffusivity_lu > 0.0);
    // The second active layer has a negative momentum fit and positive heat
    // fit. The adapter must still update its heat relaxation independently.
    for active in 20..40 {
        let cell = active + 20;
        assert_eq!(work.filter.viscosities()[active], 0.0);
        assert!(work.scalar.diffusivities()[active] > 0.0);
        assert_eq!(work.momentum[cell].to_bits(), momentum[cell].to_bits());
        assert!(work.thermal[cell] > thermal[cell]);
    }
    assert_eq!(work.momentum[..20], momentum[..20]);
    assert_eq!(work.thermal[n - 20..], thermal[n - 20..]);
    let mut heat = vec![0.0001; n];
    heat[..20].fill(0.0);
    heat[n - 20..].fill(0.0);
    let force = vec![[-2e-5, 1e-5, 0.0]; n];
    let expected = oracle
        .step_with_sources_and_transport(&heat, &force, work.relaxation(), limits(n))
        .unwrap();
    let result = state
        .step_with_sources_and_channel_sgs(
            &heat,
            &force,
            config(),
            Some(molecular),
            &mut work,
            limits(n),
            &gate,
        )
        .unwrap();
    assert_eq!(result.sgs, prepared);
    assert_eq!(result.heat, expected);
    assert_eq!(bits(&state), bits(&oracle));
    assert!(result.heat.defect.abs() < 2e-12);
    assert!((result.heat.source - 0.01).abs() < 2e-17);
    assert_eq!(state.grid.tau, original_rates);
    assert_eq!(state.tau_g, 0.9);
}

#[test]
fn g0_channel_fits_use_guo_centered_velocity_and_same_wall_temperatures() {
    let gate = CancelGate::new_clock_free();
    let mut forced = ThermalLbm3::slab(5, 4, 5, 0.8, 0.9, 0.0);
    seed(&mut forced);
    let mut unforced = forced.snapshot().0;
    for cell in 20..forced.gpop.len() - 20 {
        let u = unforced.grid.moments3(cell).1;
        let previous = unforced.grid.fext[cell];
        crate::thermal::forcing::reframe3(&mut unforced.grid.f[cell], u, previous, [0.0; 3]);
        unforced.grid.fext[cell] = [0.0; 3];
    }
    let mut a = forced.thermal_channel_sgs_workspace(usize::MAX).unwrap();
    let mut b = unforced.thermal_channel_sgs_workspace(usize::MAX).unwrap();
    forced
        .prepare_thermal_channel_sgs(config(), None, &mut a, &gate)
        .unwrap();
    unforced
        .prepare_thermal_channel_sgs(config(), None, &mut b, &gate)
        .unwrap();
    for (&actual, &expected) in a
        .filter
        .viscosities()
        .iter()
        .zip(b.filter.viscosities())
        .chain(
            a.scalar
                .diffusivities()
                .iter()
                .zip(b.scalar.diffusivities()),
        )
    {
        assert!((actual - expected).abs() < 1e-14);
    }
    // Changing a physical wall temperature must change the scalar fit while
    // leaving the same-state momentum fit independent of that scalar boundary.
    unforced.t_bottom = 0.4;
    unforced
        .prepare_thermal_channel_sgs(config(), None, &mut b, &gate)
        .unwrap();
    assert!(
        a.scalar
            .diffusivities()
            .iter()
            .zip(b.scalar.diffusivities())
            .any(|(a, b)| (a - b).abs() > 1e-8)
    );
}

#[test]
fn g0_channel_budget_heat_cap_source_cancellation_and_retry_are_atomic() {
    let gate = CancelGate::new_clock_free();
    let mut state = ThermalLbm3::slab(5, 4, 5, 0.8, 0.9, 0.0);
    seed(&mut state);
    let n = state.gpop.len();
    let before = bits(&state);
    let bytes = state
        .required_thermal_channel_sgs_workspace_bytes()
        .unwrap();
    assert!(matches!(
        state.thermal_channel_sgs_workspace(bytes - 1),
        Err(ThermalSgsError::Dynamic(
            DynamicSgsError::WorkspaceBytes { .. }
        ))
    ));
    let mut work = state.thermal_channel_sgs_workspace(bytes).unwrap();
    assert_eq!(work.workspace_bytes(), bytes);
    let mut heat = vec![0.0; n];
    let force = vec![[0.0; 3]; n];
    let mut capped = config();
    capped.heat.maximum_diffusivity_to_molecular = 0.0;
    assert!(matches!(
        state.step_with_sources_and_channel_sgs(
            &heat,
            &force,
            capped,
            None,
            &mut work,
            limits(n),
            &gate,
        ),
        Err(ThermalSgsError::DiffusivityRatio { .. })
    ));
    assert_eq!(bits(&state), before);
    heat[30] = 100.0;
    assert!(matches!(
        state.step_with_sources_and_channel_sgs(
            &heat,
            &force,
            config(),
            None,
            &mut work,
            limits(n),
            &gate,
        ),
        Err(ThermalSgsError::Source(
            ThermalSourceError::Temperature { .. }
        ))
    ));
    assert_eq!(bits(&state), before);
    heat[30] = 0.0;
    gate.request();
    assert_eq!(
        state.prepare_thermal_channel_sgs(config(), None, &mut work, &gate),
        Err(ThermalSgsError::Cancelled)
    );
    assert_eq!(bits(&state), before);
    let retry_gate = CancelGate::new_clock_free();
    let mut reference = state.snapshot().0;
    let mut fresh = reference.thermal_channel_sgs_workspace(bytes).unwrap();
    let a = state
        .step_with_sources_and_channel_sgs(
            &heat,
            &force,
            config(),
            None,
            &mut work,
            limits(n),
            &retry_gate,
        )
        .unwrap();
    let b = reference
        .step_with_sources_and_channel_sgs(
            &heat,
            &force,
            config(),
            None,
            &mut fresh,
            limits(n),
            &retry_gate,
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(bits(&state), bits(&reference));
}

#[test]
fn g1_channel_post_step_mach_refusal_keeps_original_state() {
    let gate = CancelGate::new_clock_free();
    let mut state = ThermalLbm3::slab(3, 3, 4, 0.8, 0.9, 0.0);
    let before = bits(&state);
    let n = state.gpop.len();
    let mut work = state.thermal_channel_sgs_workspace(usize::MAX).unwrap();
    let mut settings = config();
    settings.maximum_mach = 0.001;
    assert!(matches!(
        state.step_with_sources_and_channel_sgs(
            &vec![0.0; n],
            &vec![[0.002, 0.0, 0.0]; n],
            settings,
            None,
            &mut work,
            limits(n),
            &gate,
        ),
        Err(ThermalSgsError::Mach { .. })
    ));
    assert_eq!(bits(&state), before);
}

#[test]
fn g1_material_channel_adapter_conserves_si_source_and_wall_heat_in_both_orientations() {
    use crate::material::{MaterialThermalSlab3, plan_scaling_from_material};
    use fs_material::fluid_transport::FluidTransportState;
    let gate = CancelGate::new_clock_free();
    let material = FluidTransportState::try_new(2.0, 3.0, 0.02, 0.01).unwrap();
    let scale = plan_scaling_from_material(&material, 0.01, 1e-4, 0.1).unwrap();
    for (bottom, top) in [(310.0, 300.0), (300.0, 310.0)] {
        let mut slab = scale
            .thermal_slab3(MaterialThermalSlab3 {
                nx: 5,
                ny: 4,
                nz: 5,
                bottom_temperature_k: bottom,
                top_temperature_k: top,
                gravity_m_s2: 0.0,
                expansion_per_k: 0.001,
                max_cells: 140,
            })
            .unwrap();
        seed(&mut slab.lattice);
        let n = slab.lattice.gpop.len();
        let mut work = slab
            .lattice
            .thermal_channel_sgs_workspace(usize::MAX)
            .unwrap();
        let mut power = vec![100.0; n];
        power[..20].fill(0.0);
        power[n - 20..].fill(0.0);
        let result = slab
            .step_with_power_density_and_channel_sgs(
                &power,
                &vec![[0.1, 0.0, 0.0]; n],
                crate::thermal::sources::HeatSourceLimits {
                    max_cells: n,
                    temperature_range_k: [290.0, 320.0],
                    balance_tolerance_j: 1e-12,
                },
                config(),
                &mut work,
                &gate,
            )
            .unwrap();
        let prescribed = 100.0 * 100.0 * 0.01_f64.powi(3) * 1e-4;
        assert!((result.heat.source_j - prescribed).abs() < 1e-20);
        assert!(result.heat.defect_j.abs() < 1e-12);
        assert!(result.sgs.transport.momentum.maximum_eddy_viscosity_lu > 0.0);
        assert!(result.sgs.transport.maximum_eddy_diffusivity_lu > 0.0);
    }
}
