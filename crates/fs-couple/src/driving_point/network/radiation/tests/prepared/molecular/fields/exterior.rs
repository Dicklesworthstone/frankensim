//! Original complex-K exterior reaction and observations of loaded duct fields.

use super::*;
use fs_bem::helmholtz::laplace::complex::ComplexOptions;
use fs_bem::helmholtz::molecular::MolecularExterior;
use fs_bem::observation_transfers::{
    MolecularObservationFormulation, MolecularResponseLimits,
};
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchDetrend, WelchWindow};

const OMEGA: f64 = core::f64::consts::TAU * 1200.0;
const NETWORK: NetworkCoherenceBudget = NetworkCoherenceBudget {
    maximum_rank: 8,
    maximum_complex_products: 262_144,
};
const WELCH: RadiatingWelchLimits = RadiatingWelchLimits {
    maximum_segments: 8,
    maximum_mapping_products: 16_384,
    maximum_factor_bytes: 32_768,
};

fn exterior_material(source: &RelaxationState) -> MolecularExteriorMaterial<'_> {
    MolecularExteriorMaterial {
        source,
        band_rad_s: [1.0, 11_000.0],
        formulation: MolecularObservationFormulation::CentroidBurtonMiller,
        assembly: ComplexOptions {
            max_panels: 128,
            max_bytes: 1_000_000,
            max_work: 10_000_000,
        },
        response: MolecularResponseLimits {
            maximum_retained_entries: 4096,
            maximum_panel_products: 100_000,
        },
    }
}

// Real recorded pressure-generator channels, without a surrogate DFT owner.
fn recording(drives: &[[C64; 2]]) -> WelchCrossSpectrum {
    const N: usize = 32;
    let rows: Vec<_> = drives.iter().flat_map(|drive| {
        (0..N).map(move |j| {
            let phase = core::f64::consts::TAU * 3.0 * j as f64 / N as f64;
            drive.map(|z| z.re * phase.cos() + z.im * phase.sin())
        })
    }).collect();
    WelchCrossSpectrum::estimate_interleaved(
        &rows,
        1.0 / 12_800.0,
        WelchConfig {
            segment_samples: N,
            hop_samples: N,
            window: WelchWindow::Rectangular,
            detrend: WelchDetrend::None,
        },
        WelchBudget {
            maximum_channels: 2,
            maximum_samples_per_channel: rows.len(),
            maximum_transforms: 2 * drives.len(),
            maximum_workspace_bytes: 1_000_000,
        },
        || false,
    ).unwrap()
}

#[test]
fn original_molecular_exterior_closes_the_same_loaded_material_wall_fields() {
    let model = regional();
    let gate = CancelGate::new_clock_free();
    let wall_source = wall();
    let bindings = walls(&model, &wall_source);
    let original_wall = bindings[0].1.clone();
    let edge = bindings[0].0.clone();
    let molecular = model.bind_molecular_with_compliant_walls(
        edges(&model, 2, 4.0), bindings, 64, &gate,
    ).unwrap();
    let atmosphere = state(model.gas(), 2, 8.0);
    let material = exterior_material(&atmosphere);
    let (probes, cells) = probes(&model);
    let bank = bank(OMEGA);
    let prepared = molecular.prepare_fields_with_exterior_material(
        &request(OMEGA, &bank, &probes), material, &gate,
    ).unwrap();
    let radiating = prepared.radiating();
    assert!(core::ptr::eq(radiating.exterior_material().unwrap().source, &atmosphere));
    assert!(core::ptr::eq(
        prepared.interior().network().edge(&edge).unwrap().compliant_wall().unwrap(),
        original_wall.as_ref(),
    ));
    assert!(radiating.far_field().is_none());
    assert!(radiating.boundary().z[1].abs() > 0.0);
    let calls = wall_source.calls.load(Ordering::Relaxed);
    assert_eq!(calls, model.driving_point().network().edges.len());

    let amplitudes = radiating.default_excitations();
    let mut actual = radiating.response_buffer();
    radiating.drive_into(amplitudes, &mut actual, &gate).unwrap();
    let original = prepared.interior().fields().port_response();
    let mut local = original.response_buffer();
    original.drive_into(amplitudes, &mut local).unwrap();
    assert_eq!(actual.network().unwrap().fields.port_pressures, local.fields.port_pressures);
    assert_eq!(actual.network().unwrap().fields.port_inflows, local.fields.port_inflows);

    // Independently assemble the complete physical panel velocity, then solve
    // it directly. This does not use the prepared unit-opening transfer bank.
    let outward = actual.outward_flows_m3_s().unwrap();
    let mut velocity = vec![C64::ZERO; model.surface().areas().len()];
    for ((patch, &area), &flow) in model.patches().iter()
        .zip(model.opening_areas_m2()).zip(outward)
    {
        for &panel in &patch.panels {
            velocity[panel] = flow.scale(area.recip());
        }
    }
    let direct = MolecularExterior::new(atmosphere.clone(), material.band_rad_s).unwrap();
    let solved = direct.solve_batch(
        model.surface(), OMEGA, &[velocity.as_slice()], Formulation::BurtonMiller,
        &material.assembly, || false,
    ).unwrap();
    let fields = solved[0].observe(
        model.surface(), model.receivers(), model.clearance_m, || false,
    ).unwrap();
    let euler = C64::new(0.0, -1.0 / (OMEGA * atmosphere.gas().density));
    for (a, b) in actual.receiver_fields().unwrap().iter().zip(fields) {
        close(a.pressure_pa, b.pressure_pa);
        for axis in 0..3 {
            close(a.particle_velocity_m_s[axis], euler * b.gradient_pa_per_m[axis]);
        }
    }
    for patch in model.patches() {
        let area: f64 = patch.panels.iter().map(|&i| model.surface().areas()[i]).sum();
        let pressure = patch.panels.iter().fold(C64::ZERO, |sum, &i| {
            sum + solved[0].pressure_pa()[i].scale(model.surface().areas()[i] / area)
        });
        let index = radiating.ports().iter().position(|label| label == &patch.port).unwrap();
        close(pressure, local.fields.port_pressures[index]);
    }
    real_close(actual.exterior_surface_power_w().unwrap(), solved[0].boundary_work_w());
    let exterior_work = radiating.molecular_exterior().unwrap().power_w(outward, || false).unwrap();
    real_close(exterior_work.driving_power_w, solved[0].boundary_work_w());
    real_close(exterior_work.exterior_power_w, exterior_work.driving_power_w);
    assert_eq!(exterior_work.absorbed_wall_power_w, 0.0);
    let absorption = prepared.interior().prepare_cylinder_absorption(
        &cells, cells.len(), || false,
    ).unwrap();
    let losses = absorption.sample(amplitudes, || false).unwrap();
    let internal: f64 = losses.iter().map(|loss| loss.total_removed_power_w()).sum();
    assert!(losses.iter().any(|loss| loss.wall_work_w > 0.0));
    real_close(internal, local.fields.network_power_w);
    real_close(local.generator_power_w,
        local.source_dissipation_w + internal + exterior_work.exterior_power_w);
    assert_eq!(wall_source.calls.load(Ordering::Relaxed), calls);

    // Only the external reservoir rates change. The load and solved mouth
    // flows must change before observation, not merely receive an output filter.
    let fast = state(model.gas(), 2, 0.05);
    let changed = molecular.prepare_fields_with_exterior_material(
        &request(OMEGA, &bank, &probes), exterior_material(&fast), &gate,
    ).unwrap();
    let mut other = changed.radiating().response_buffer();
    changed.radiating().drive_into(amplitudes, &mut other, &gate).unwrap();
    assert_ne!(radiating.boundary().source_id, changed.radiating().boundary().source_id);
    assert!((radiating.boundary().z[0] - changed.radiating().boundary().z[0]).abs()
        > 1e-6 * radiating.boundary().z[0].abs());
    assert!(outward.iter().zip(other.outward_flows_m3_s().unwrap())
        .any(|(&a, &b)| (a - b).abs() > 1e-6 * a.abs()));
}

#[test]
fn empty_reservoir_limit_and_recorded_coherence_reuse_the_native_observers() {
    let model = regional();
    let gate = CancelGate::new_clock_free();
    let molecular = model.bind_molecular(edges(&model, 2, 4.0), 64, &gate).unwrap();
    let (probes, _) = probes(&model);
    let bank = bank(OMEGA);
    let empty = state(model.gas(), 0, 1.0);
    let material = molecular.prepare_fields_with_exterior_material(
        &request(OMEGA, &bank, &probes), exterior_material(&empty), &gate,
    ).unwrap();
    let ordinary = molecular.prepare_fields(&request(OMEGA, &bank, &probes), &gate).unwrap();
    for (&a, &b) in material.radiating().boundary().z.iter().zip(&ordinary.radiating().boundary().z) {
        close(a, b);
    }
    let mut a = material.radiating().response_buffer();
    let mut b = ordinary.radiating().response_buffer();
    material.radiating().drive_into(material.radiating().default_excitations(), &mut a, &gate).unwrap();
    ordinary.radiating().drive_into(ordinary.radiating().default_excitations(), &mut b, &gate).unwrap();
    for (a, b) in a.receiver_fields().unwrap().iter().zip(b.receiver_fields().unwrap()) {
        close(a.pressure_pa, b.pressure_pa);
        for axis in 0..3 { close(a.particle_velocity_m_s[axis], b.particle_velocity_m_s[axis]); }
    }
    assert!(ordinary.radiating().exterior_material().is_none());

    let atmosphere = state(model.gas(), 2, 8.0);
    let prepared = molecular.prepare_fields_with_exterior_material(
        &request(OMEGA, &bank, &probes), exterior_material(&atmosphere), &gate,
    ).unwrap();
    let owner = prepared.radiating();
    let drives = [
        [C64::new(0.8, -0.2), C64::new(-0.3, 0.7)],
        [C64::new(0.1, 0.4), C64::new(0.2, -0.6)],
    ];
    let factors: Vec<_> = drives.iter().flatten().map(|z| z.scale(0.5_f64.sqrt())).collect();
    let spectrum = recording(&drives);
    let measured = owner.prepare_welch_bin(
        &spectrum, 3, owner.source_coordinates(), WELCH, &gate,
    ).unwrap();
    let mut actual = owner.statistics_buffer();
    measured.evaluate_into(NETWORK, FIELD_BUDGET, &mut actual, &gate).unwrap();
    let mut reference = owner.statistics_buffer();
    owner.coherence_into(&factors, 2, NETWORK, FIELD_BUDGET, &mut reference, &gate).unwrap();
    let mut independent_pressure = vec![0.0; model.receivers().len()];
    for drive in &drives {
        let mut field = owner.response_buffer();
        owner.drive_into(drive, &mut field, &gate).unwrap();
        for (power, value) in independent_pressure.iter_mut().zip(field.receiver_fields().unwrap()) {
            *power += 0.25 * value.pressure_pa.norm_sq();
        }
    }
    for ((a, b), expected) in actual.receiver_fields().unwrap().iter()
        .zip(reference.receiver_fields().unwrap()).zip(independent_pressure)
    {
        real_close(a.mean_square_pressure_pa2, expected);
        real_close(a.mean_square_pressure_pa2, b.mean_square_pressure_pa2);
        for axis in 0..3 { close(a.complex_intensity_w_m2[axis], b.complex_intensity_w_m2[axis]); }
    }
    let patterns = [
        FirstOrder::new(0.4, [1.0, 0.0, 0.0]).unwrap(),
        FirstOrder::new(0.2, [0.0, 1.0, 0.0]).unwrap(),
    ];
    let mut microphones = owner.microphones(&patterns).unwrap();
    let expected = microphones.coherence(&factors, 2, FIELD_BUDGET, &gate).unwrap().to_vec();
    for (&a, b) in microphones.welch_bin(&measured, FIELD_BUDGET, &gate).unwrap().iter().zip(expected) {
        close(a, b);
    }
    let work = actual.network().unwrap().work().unwrap();
    real_close(work.generator_power_w,
        work.source_dissipation_w + work.network_power_w + work.external_load_power_w);
    real_close(actual.exterior_surface_power_w().unwrap(), work.external_load_power_w);
    assert!(owner.far_field().is_none());
}

#[test]
fn exterior_domains_and_stationary_band_source_identity_are_binding() {
    let model = regional();
    let gate = CancelGate::new_clock_free();
    let molecular = model.bind_molecular(edges(&model, 2, 4.0), 64, &gate).unwrap();
    let atmosphere = state(model.gas(), 2, 8.0);
    let other_source = atmosphere.clone();
    let hot = GasState::try_new(
        &fs_material::gas::GasSpec::dry_air_ussa1976(), 310.0, model.gas().pressure,
    ).unwrap();
    let foreign_gas = state(hot, 2, 8.0);
    let (probes, _) = probes(&model);
    let frequencies = [OMEGA, core::f64::consts::TAU * 1600.0];
    let banks = frequencies.map(bank);
    let request0 = request(frequencies[0], &banks[0], &probes);
    assert!(matches!(molecular.prepare_fields_with_exterior_material(
        &request0, exterior_material(&foreign_gas), &gate,
    ), Err(ExteriorNetworkError::Invalid(_))));
    for case in 0..3 {
        let mut material = exterior_material(&atmosphere);
        match case {
            0 => material.band_rad_s = [OMEGA * 1.1, OMEGA * 1.2],
            1 => material.assembly.max_bytes = 0,
            _ => material.response.maximum_retained_entries = 0,
        }
        assert!(matches!(molecular.prepare_fields_with_exterior_material(
            &request0, material, &gate,
        ), Err(ExteriorNetworkError::MolecularExterior(_))));
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(molecular.prepare_fields_with_exterior_material(
        &request0, exterior_material(&atmosphere), &cancelled,
    ), Err(ExteriorNetworkError::Cancelled)));

    let owners: Vec<_> = frequencies.iter().zip(&banks).map(|(&omega, bank)| {
        molecular.prepare_fields_with_exterior_material(
            &request(omega, bank, &probes), exterior_material(&atmosphere), &gate,
        ).unwrap()
    }).collect();
    let foreign = molecular.prepare_fields_with_exterior_material(
        &request(frequencies[1], &banks[1], &probes), exterior_material(&other_source), &gate,
    ).unwrap();
    let ordinary = molecular.prepare_fields(
        &request(frequencies[1], &banks[1], &probes), &gate,
    ).unwrap();
    let spectrum = recording(&[[C64::ONE, C64::new(0.2, 0.4)]]);
    let limits = RadiatingWelchBandLimits {
        maximum_bins: 2,
        maximum_total_factor_bytes: 100_000,
        maximum_total_mapping_products: 100_000,
    };
    let bins = [
        owners[0].radiating().prepare_welch_bin(
            &spectrum, 3, owners[0].radiating().source_coordinates(), WELCH, &gate,
        ).unwrap(),
        owners[1].radiating().prepare_welch_bin(
            &spectrum, 4, owners[1].radiating().source_coordinates(), WELCH, &gate,
        ).unwrap(),
    ];
    assert!(PreparedRadiatingWelchBand::stationary(&bins, limits, &gate).is_ok());
    // The second recorded bin is silent, but source identity still matters.
    // Numerical equality of a cloned law cannot authorize a different owner.
    for replacement in [foreign.radiating(), ordinary.radiating()] {
        let mixed = [
            owners[0].radiating().prepare_welch_bin(
                &spectrum, 3, owners[0].radiating().source_coordinates(), WELCH, &gate,
            ).unwrap(),
            replacement.prepare_welch_bin(
                &spectrum, 4, replacement.source_coordinates(), WELCH, &gate,
            ).unwrap(),
        ];
        assert!(matches!(PreparedRadiatingWelchBand::stationary(&mixed, limits, &gate),
            Err(ExteriorNetworkError::Invalid(_))));
    }
}
