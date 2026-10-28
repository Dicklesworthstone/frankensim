//! Actual detached-interface heat, sourced surface stress and native evolution.
//! Authored liquid/gas data verify the coupling, not experimental calibration.

use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceConductionControls, FreeSurfaceInterfaceContact, FreeSurfaceInterfaceSolve,
    LocalSurfaceBoundaryEvolution,
};

#[path = "../../../nsf3_thermochemical/source.rs"]
mod gas_source;

#[cfg(feature = "implicit-enthalpy-coupling")]
#[path = "surface_heat/implicit.rs"]
mod implicit;

fn heat_controls() -> FreeSurfaceConductionControls {
    FreeSurfaceConductionControls {
        maximum_faces: 12000,
        maximum_substeps: 10000,
        maximum_dt_s: 0.00025,
        energy_tolerance_j: 1e-7,
    }
}

fn uniform_solve() -> FreeSurfaceInterfaceSolve {
    FreeSurfaceInterfaceSolve {
        local_boundary_evolution: LocalSurfaceBoundaryEvolution::UniformReservoir,
        ..FreeSurfaceInterfaceSolve::default()
    }
}

#[test]
fn g1_local_solved_heat_uses_live_patches_shared_sigma_and_one_heat_account() {
    let slope = 1e-6;
    let gas = gas_source::model(0.0, 250000.0, 1000.0);
    let point = gas_source::point(450.0, 100000.0);
    let controls = FreeSurfaceConductionControls {
        // The actual thermal CFL still bounds this interval. This regression
        // targets moving geometry; substep accumulation has separate tests.
        maximum_dt_s: DT,
        ..heat_controls()
    };
    // A small prescribed translation exposes new patch owners within the two
    // actual intervals; capillary startup alone barely moves this large drop.
    with_moving_drop(slope, [0.01, 0.0, 0.0], [1.0, 0.0, 0.0], |thermal, gate| {
        let mut changed_owners = false;
        for _ in 0..2 {
            let patches = thermal.local_interface_patches(gate).unwrap();
            let owners: std::collections::BTreeSet<_> = patches.iter().map(|p| p.liquid).collect();
            let boundaries: Vec<_> = patches
                .iter()
                .map(|patch| FreeSurfaceInterfaceContact {
                    surface: patch.with_exchange([0.25 * DX, 0.5 * DX], 0.0, 400.0),
                    model: &gas,
                    context: &point,
                })
                .collect();
            let mass = thermal.flow().ledger_mass_kg().unwrap();
            let energy = thermal.total_enthalpy_j().unwrap();
            let report = thermal
                .step_with_surface_exchange(
                    &[],
                    &[],
                    controls,
                    &mut [],
                    &[],
                    &[],
                    &boundaries,
                    uniform_solve(),
                    false,
                    gate,
                )
                .unwrap()
                .conduction;
            let exchange = report.surface.as_ref().unwrap();
            assert!(exchange.heat_from_gas_j > 0.0);
            assert_eq!(exchange.heat_from_radiation_j, 0.0);
            near(report.flow.external_heat_j, exchange.heat_from_gas_j, 2e-9);
            assert!(
                (thermal.total_enthalpy_j().unwrap() - energy - exchange.heat_from_gas_j).abs()
                    < 1e-7
            );
            near(thermal.flow().ledger_mass_kg().unwrap(), mass, 3e-12);

            let live = thermal.local_interface_patches(gate).unwrap();
            assert_eq!(exchange.native_endpoint.len(), live.len());
            let new_owners: std::collections::BTreeSet<_> = live.iter().map(|p| p.liquid).collect();
            changed_owners |= owners != new_owners;
            let samples: std::collections::BTreeMap<_, _> = exchange
                .native_endpoint
                .iter()
                .map(|s| (s.surface.face.liquid, s))
                .collect();
            let mut net_x = 0.0;
            for patch in &live {
                let sample = samples[&patch.liquid];
                assert_eq!(sample.surface.face.gas, patch.liquid);
                // Both APIs convert the shared lattice patch to SI. Compare
                // their roundoff at the cell-area scale, including tiny cuts;
                // a source-frame change must preserve the actual geometry.
                assert!(
                    (sample.surface.face.area_m2 - patch.area_m2).abs()
                        <= 64.0 * f64::EPSILON * DX * DX
                );
                assert!(patch.liquid_boundary_area_m2 >= patch.area_m2);
                // Independent two-layer solution. The original authored liquid
                // has k=0.5 W/(m K); gas k depends on pressure, not film T here.
                let kl = 0.5;
                let kg = 0.026 * 1000.0 * (0.5 + 99000.0 / 249000.0);
                let gl = kl / (0.25 * DX);
                let gg = kg / (0.5 * DX);
                let ts = (gl * sample.liquid_temperature_k + gg * 450.0) / (gl + gg);
                near(sample.surface_temperature_k, ts, 3e-11);
                let sigma = 0.001 + slope * (ts - 400.0);
                near(
                    sample.surface_tension.as_ref().unwrap().value_si(),
                    sigma,
                    3e-12,
                );
                let [x, y, z] = patch.liquid;
                near(
                    thermal
                        .flow()
                        .lattice()
                        .surface_tension_at(x, y, z)
                        .unwrap(),
                    sigma * DT * DT / (1000.0 * DX.powi(3)),
                    3e-12,
                );
                let bulk_sigma = thermal
                    .material_at(patch.liquid)
                    .unwrap()
                    .material()
                    .property(SURFACE_TENSION_PROPERTY)
                    .unwrap()
                    .value_si();
                assert!((sigma - bulk_sigma).abs() > 1e-6);
                let force = thermal.flow().lattice().force_density_at(x, y, z).unwrap();
                let fill = thermal
                    .flow()
                    .lattice()
                    .sample_at(x, y, z)
                    .unwrap()
                    .fill_fraction;
                assert!(
                    (0..3)
                        .map(|d| force[d] * patch.outward_normal[d])
                        .sum::<f64>()
                        .abs()
                        < 2e-14
                );
                net_x += fill * force[0];
            }
            assert!(
                net_x > 0.0,
                "the solved surface temperature gradient must drive native tangential momentum"
            );
        }
        assert!(
            changed_owners,
            "this evolving heat regression must exercise newly exposed local patches"
        );
    });
}

#[test]
fn g4_constant_sigma_local_heat_owns_live_geometry_without_locking_plain_steps() {
    // Zero heat and a weak bath conductance keep this a single thermal
    // interval. The actual translation still creates/removes native patches.
    let gas = gas_source::model(0.0, 250000.0, 1e-6);
    let point = gas_source::point(400.0, 100000.0);
    with_moving_drop(0.0, [0.01, 0.0, 0.0], [0.0; 3], |thermal, gate| {
        let patches = thermal.local_interface_patches(gate).unwrap();
        let owners: std::collections::BTreeSet<_> = patches.iter().map(|p| p.liquid).collect();
        let boundaries: Vec<_> = patches
            .iter()
            .map(|patch| FreeSurfaceInterfaceContact {
                surface: patch.with_exchange([0.25 * DX, 0.5 * DX], 0.0, 400.0),
                model: &gas,
                context: &point,
            })
            .collect();
        let controls = FreeSurfaceConductionControls {
            maximum_substeps: 1,
            maximum_dt_s: DT,
            ..heat_controls()
        };
        let populations = thermal.flow().lattice().grid().state_bits();
        let energy = thermal.total_enthalpy_j().unwrap();
        let mass = thermal.flow().ledger_mass_kg().unwrap();
        let error = thermal
            .step_with_surface_exchange(
                &[],
                &[],
                controls,
                &mut [],
                &[],
                &[],
                &boundaries,
                FreeSurfaceInterfaceSolve::default(),
                false,
                gate,
            )
            .unwrap_err();
        assert!(matches!(
            error,
            FreeSurfaceEnthalpyError::Invalid("fixed local heat owners changed during native flow")
        ));
        assert_eq!(populations, thermal.flow().lattice().grid().state_bits());
        assert_eq!(
            energy.to_bits(),
            thermal.total_enthalpy_j().unwrap().to_bits()
        );
        assert_eq!(
            mass.to_bits(),
            thermal.flow().ledger_mass_kg().unwrap().to_bits()
        );

        let report = thermal
            .step_with_surface_exchange(
                &[],
                &[],
                controls,
                &mut [],
                &[],
                &[],
                &boundaries,
                uniform_solve(),
                false,
                gate,
            )
            .unwrap()
            .conduction;
        let surface = report.surface.as_ref().unwrap();
        let live = thermal.local_interface_patches(gate).unwrap();
        assert_eq!(surface.native_endpoint.len(), live.len());
        assert_ne!(
            owners,
            live.iter()
                .map(|p| p.liquid)
                .collect::<std::collections::BTreeSet<_>>()
        );
        assert_eq!(surface.heat_from_gas_j, 0.0);
        assert_eq!(surface.heat_from_radiation_j, 0.0);
        for sample in &surface.native_endpoint {
            assert_eq!(sample.surface_tension.as_ref().unwrap().value_si(), 0.001);
        }
        assert!((thermal.total_enthalpy_j().unwrap() - energy).abs() < 1e-7);
        // A true scalar sigma needs no interface temperature for its force;
        // requesting one diagnostic heat step must not permanently bind heat.
        let next = thermal.step(gate).unwrap();
        assert_eq!(next.external_heat_j, 0.0);
        near(next.mass_after_kg, mass, 3e-12);
        assert!((next.enthalpy_after_j - energy).abs() < 1e-7);
    });
}

#[test]
fn g4_local_solved_heat_rejects_incomplete_geometry_and_nonuniform_continuation_atomically() {
    let gas = gas_source::model(0.0, 250000.0, 1000.0);
    let point = gas_source::point(450.0, 100000.0);
    let different = gas_source::point(449.0, 100000.0);
    for slope in [0.0, 1e-6] {
        with_drop(slope, |thermal, gate| {
            let boundaries: Vec<_> = thermal
                .local_interface_patches(gate)
                .unwrap()
                .into_iter()
                .map(|patch| FreeSurfaceInterfaceContact {
                    surface: patch.with_exchange([0.25 * DX, 0.5 * DX], 0.0, 400.0),
                    model: &gas,
                    context: &point,
                })
                .collect();
            let populations = thermal.flow().lattice().grid().state_bits();
            let energy = thermal.total_enthalpy_j().unwrap().to_bits();
            let mass = thermal.flow().ledger_mass_kg().unwrap().to_bits();
            let [x, y, z] = TANGENT_POLE;
            let force = thermal.flow().lattice().force_density_at(x, y, z);
            let sigma = thermal.flow().lattice().surface_tension_at(x, y, z);
            for kind in 0..6 {
                let mut bad = boundaries.clone();
                match kind {
                    0 => {
                        bad.pop();
                    }
                    1 => {
                        bad[1] = bad[0];
                    }
                    2 => {
                        bad[0].surface.face.area_m2 *= 0.99;
                    }
                    3 => {
                        bad[0].surface.face.gas[0] += 1;
                    }
                    4 => {
                        bad[0].surface.face.distances_m[0] = 0.0;
                    }
                    _ => {
                        bad[0].context = &different;
                    }
                }
                assert!(
                    thermal
                        .step_with_surface_exchange(
                            &[],
                            &[],
                            heat_controls(),
                            &mut [],
                            &[],
                            &[],
                            &bad,
                            uniform_solve(),
                            false,
                            gate,
                        )
                        .is_err()
                );
                assert_eq!(populations, thermal.flow().lattice().grid().state_bits());
                assert_eq!(energy, thermal.total_enthalpy_j().unwrap().to_bits());
                assert_eq!(mass, thermal.flow().ledger_mass_kg().unwrap().to_bits());
                assert_eq!(force, thermal.flow().lattice().force_density_at(x, y, z));
                assert_eq!(sigma, thermal.flow().lattice().surface_tension_at(x, y, z));
            }
        });
    }
}
