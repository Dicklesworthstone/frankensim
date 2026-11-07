//! G1/G3/G4: physical heated parcels anchored at accepted cavity geometry.

use super::*;
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_solid::hyper3::pressure::cavity::{AdiabaticCavity, AdiabaticCavityGas, CavityGasEpoch};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn controls() -> GasEntropySolve {
    GasEntropySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 650.0,
        entropy_tolerance_j_kg_k: 1e-10,
        maximum_iterations: 128,
    }
}

fn displacement(nodes: &[[f64; 3]], stretch: f64) -> Vec<f64> {
    nodes.iter().flatten().map(|x| (stretch - 1.0) * x).collect()
}

/// Independent variable-Cv entropy integral; no native entropy solver is used.
fn ideal_temperature(initial: f64, volume_ratio: f64, slope: f64) -> f64 {
    let (mut lo, mut hi) = (200.0, 650.0);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        let residual = 2.5 * (mid / initial).ln() + slope * (mid - initial)
            + volume_ratio.ln();
        if residual < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn heated_deformed_epoch_starts_the_actual_isentrope_and_anchor_work() {
    with_cx(|_, cx| {
        let (nodes, tets) = shell();
        let cavity = Arc::new(cavity(&nodes, &tets, cx));
        let anchor_stretch = 0.75_f64;
        let target_stretch = 0.72_f64;
        let anchor_u = displacement(&nodes, anchor_stretch);
        let target_u = displacement(&nodes, target_stretch);
        let anchor_volume = cavity.evaluate(&anchor_u, cx).unwrap().volume_m3();
        let slope = 0.001;
        let original = gas_fixture::gas(cavity.reference_volume_m3(), 100000.0, slope, 1e24);
        let mut compressed = original.clone();
        compressed.set_volume_isentropic(anchor_volume, controls(), || false).unwrap();
        let before_heat = compressed.state().state().temperature;
        let (heated, supplied_heat) = compressed.at_temperature(500.0).unwrap();
        let capacity_scale = original.mass_kg() * R / 0.028;
        near(supplied_heat,
            capacity_scale * (500.0 - before_heat) * (2.5 + 0.5 * slope * (500.0 + before_heat)),
            2e-7);
        assert!(supplied_heat > 1000.0);
        assert_eq!(heated.volume_m3().to_bits(), anchor_volume.to_bits());
        let saved = heated.clone();
        let epoch = CavityGasEpoch::try_new(
            Arc::clone(&cavity), &anchor_u, heated, cx,
        ).unwrap();
        assert_eq!(epoch.gas(), &saved);
        assert!(Arc::ptr_eq(epoch.cavity(), &cavity));
        assert!(Arc::ptr_eq(epoch.clone().cavity(), epoch.cavity()));
        assert_eq!(epoch.anchor_volume_m3().to_bits(), anchor_volume.to_bits());
        let narrow = GasEntropySolve {
            minimum_temperature_k: 480.0,
            maximum_temperature_k: 560.0,
            ..controls()
        };
        // A fictitious return to undeformed volume is outside this physical
        // interval's temperature chart. The accepted-volume reference needs no
        // such state or extrapolation.
        let mut false_rebase = saved.clone();
        assert!(false_rebase.set_volume_isentropic(
            cavity.reference_volume_m3(), narrow, || false,
        ).is_err());
        assert_eq!(false_rebase, saved);
        let ambient = 20000.0;
        let binding = AdiabaticCavity {
            cavity: cavity.as_ref(), reference: &epoch,
            exterior_pressure_pa: ambient, controls: narrow,
        };
        let start = binding.evaluate(&anchor_u, cx).unwrap();
        assert_eq!(start.gas.gas(), epoch.gas());
        assert_eq!(start.potential_change_j, 0.0);
        // Repeated physical heat probes retain one admitted immutable anchor.
        let (hotter, extra_heat) = epoch.gas().at_temperature(520.0).unwrap();
        assert!(extra_heat > 0.0);
        let warmer = epoch.with_gas_at_anchor(hotter.clone(), cx).unwrap();
        assert!(Arc::ptr_eq(warmer.cavity(), epoch.cavity()));
        assert_eq!(warmer.anchor_volume_m3(), epoch.anchor_volume_m3());
        assert_eq!(warmer.gas(), &hotter);
        let warm_start = cavity.adiabatic_gas_load(
            &anchor_u, &warmer, ambient, narrow, cx,
        ).unwrap();
        assert_eq!(warm_start.potential_change_j, 0.0);
        assert_eq!(warm_start.gas.gas(), &hotter);
        assert_eq!(epoch.gas(), &saved);
        let load = binding.evaluate(&target_u, cx).unwrap();
        let (volume, gradient) = tetra_volume_gradient(current_inner(&nodes, &target_u));
        let temperature = ideal_temperature(
            500.0, (target_stretch / anchor_stretch).powi(3), slope,
        );
        let pressure = capacity_scale * temperature / volume;
        let stored = capacity_scale * (temperature - 500.0)
            * (2.5 + 0.5 * slope * (temperature + 500.0));
        near(load.gas.gas().state().state().temperature, temperature, 2e-8);
        near(load.net_pressure_pa, pressure - ambient, 2e-5);
        near(load.potential_change_j, stored + ambient * (volume - anchor_volume), 2e-6);
        near(start.gas.stored_energy_change_to(&load.gas, || false).unwrap(), stored, 2e-6);
        assert_eq!(load.gas.gas().mass_kg().to_bits(), original.mass_kg().to_bits());
        assert_eq!(load.gas.gas().model(), original.model());
        assert_eq!(load.gas.gas().point().axes().get("material-lot"), Some(&47.0));
        let gamma = (3.5 + slope * temperature) / (2.5 + slope * temperature);
        for local in 0..12 {
            near(load.force_n[12 + local], (pressure - ambient) * gradient[local], 2e-5);
            near(load.force_tangent.acoustic_rank_one_vector()[12 + local],
                (gamma * pressure / volume).sqrt() * gradient[local], 2e-7);
        }
        // One directional derivative checks both the acoustic and geometric
        // force terms, plus the new accepted-anchor potential's work sign.
        let direction: Vec<_> = (0..24).map(|i| ((i * 5 % 13) as f64 - 6.0) / 30.0).collect();
        let h = 2e-5;
        let plus: Vec<_> = target_u.iter().zip(&direction).map(|(u, d)| u + h * d).collect();
        let minus: Vec<_> = target_u.iter().zip(&direction).map(|(u, d)| u - h * d).collect();
        let plus = binding.evaluate(&plus, cx).unwrap();
        let minus = binding.evaluate(&minus, cx).unwrap();
        let work = pairing(&load.force_n, &direction);
        near((plus.potential_change_j - minus.potential_change_j) / (2.0 * h),
            -work, 2e-6 * work.abs().max(1.0));
        let tangent = load.force_tangent.apply(&direction, cx).unwrap();
        let scale = tangent.iter().map(|v| v.abs()).fold(1.0_f64, f64::max);
        for (i, value) in tangent.iter().enumerate() {
            near(*value, (plus.force_n[i] - minus.force_n[i]) / (2.0 * h), 2e-6 * scale);
        }
        assert_eq!(epoch.gas(), &saved);
        assert_eq!(compressed.mass_kg().to_bits(), original.mass_kg().to_bits());
    });
}

#[test]
fn epochs_require_explicit_geometry_rebinding_and_refuse_atomically() {
    with_cx(|gate, cx| {
        let (nodes, tets) = shell();
        let cavity = Arc::new(cavity(&nodes, &tets, cx));
        let anchor_u = displacement(&nodes, 0.8);
        let after_u = displacement(&nodes, 0.78);
        let volume = cavity.evaluate(&anchor_u, cx).unwrap().volume_m3();
        let gas = gas_fixture::gas(volume, 240000.0, 0.0, 0.0);
        let saved = gas.clone();
        let epoch = CavityGasEpoch::try_new(
            Arc::clone(&cavity), &anchor_u, gas.clone(), cx,
        ).unwrap();
        let start = cavity.adiabatic_gas_load(&anchor_u, &epoch, 10000.0, controls(), cx).unwrap();
        let after = cavity.adiabatic_gas_load(&after_u, &epoch, 10000.0, controls(), cx).unwrap();
        // Ordinary parcels still require the undeformed volume, and even an
        // equal-volume output epoch cannot silently become a new reference.
        assert!(cavity.adiabatic_gas_load(&anchor_u, &gas, 0.0, controls(), cx).is_err());
        for endpoint in [&start.gas, &after.gas] {
            assert!(cavity.adiabatic_gas_load(
                &after_u, endpoint, 0.0, controls(), cx,
            ).is_err());
            assert!(endpoint.with_gas_at_anchor(saved.clone(), cx).is_err());
        }
        assert!(epoch.with_gas_at_anchor(after.gas.gas().clone(), cx).is_err());
        let impostor = Arc::new(cavity.as_ref().clone());
        assert!(impostor.adiabatic_gas_load(
            &after_u, &epoch, 0.0, controls(), cx,
        ).is_err());
        // Explicitly rebinding an unchanged physical endpoint is legal, and
        // resets only its mechanical potential anchor.
        let rebound = CavityGasEpoch::try_new(
            Arc::clone(&cavity), &after_u, after.gas.clone().into_gas(), cx,
        ).unwrap();
        let rebound_start = cavity.adiabatic_gas_load(
            &after_u, &rebound, 10000.0, controls(), cx,
        ).unwrap();
        assert_eq!(rebound_start.gas.gas(), after.gas.gas());
        assert_eq!(rebound_start.potential_change_j, 0.0);
        assert_eq!(rebound_start.force_n, after.force_n);
        assert!(CavityGasEpoch::try_new(
            Arc::clone(&cavity), &after_u, gas.clone(), cx,
        ).is_err());
        assert!(CavityGasEpoch::try_new(
            Arc::clone(&cavity), &[0.0; 23], gas.clone(), cx,
        ).is_err());
        // The anchor admits the entire solid, including outer elements away
        // from the selected gas wall.
        let mut inverted = anchor_u.clone();
        inverted[3] = -4.0;
        assert!(CavityGasEpoch::try_new(
            Arc::clone(&cavity), &inverted, gas.clone(), cx,
        ).is_err());
        let invalid_controls = GasEntropySolve {
            maximum_temperature_k: 301.0, ..controls()
        };
        assert!(cavity.adiabatic_gas_load(
            &after_u, &epoch, 10000.0, invalid_controls, cx,
        ).is_err());
        assert_eq!(epoch.gas(), &saved);
        gate.request();
        assert!(matches!(epoch.with_gas_at_anchor(gas.clone(), cx),
            Err(HyperTetError::Cancelled)));
        assert!(matches!(CavityGasEpoch::try_new(
            Arc::clone(&cavity), &anchor_u, gas.clone(), cx,
        ), Err(HyperTetError::Cancelled)));
        assert!(matches!(cavity.adiabatic_gas_load(
            &after_u, &epoch, 0.0, controls(), cx,
        ), Err(HyperTetError::Cancelled)));
        assert_eq!(epoch.gas(), &saved);
        assert_eq!(gas, saved);
    });
}
