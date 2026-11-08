use super::*;
use crate::chemistry::diffuse_recombination::DiffuseRecombinationBeam;
use crate::chemistry::ionization_network::transport::*;

fn beam(energy_ev: f64, flux: f64, mu: f64) -> DiffuseRecombinationBeam {
    DiffuseRecombinationBeam {
        photons: PhotonNumberFluxV1 {
            photon_energy: Energy::new(energy_ev * EV),
            number_flux_per_m2_s: flux,
        },
        direction_cosine: mu,
    }
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    let scale = actual.abs().max(expected.abs()).max(f64::MIN_POSITIVE);
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance * scale,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn scaled(state: &FrozenIonizedGasStateV1, fraction: f64) -> FrozenIonizedGasStateV1 {
    FrozenIonizedGasStateV1::new(
        state.thermodynamics(),
        &state
            .amounts()
            .iter()
            .map(|n| Amount::new(n.value() * fraction))
            .collect::<Vec<_>>(),
        Volume::new(state.volume().value() * fraction),
        state.heavy_temperature(),
        state.electron_temperature(),
        1e-12,
        || false,
    )
    .unwrap()
}

// Independent constant-cross-section optical depths from the manufactured
// fixture. The oracle calls neither production transfer nor cross-section code.
fn channel_depths(state: &FrozenIonizedGasStateV1, area: f64, energy_ev: f64) -> [f64; 3] {
    core::array::from_fn(|edge| {
        if energy_ev < GAPS_EV[edge] {
            0.0
        } else {
            SIGMA[edge] * NA * state.amounts()[EDGES[edge][0]].value() / area
        }
    })
}

#[test]
fn shared_photons_follow_summed_opacity_from_thin_to_opaque() {
    let (state, sources) = fixture(false, [0.0; 3]);
    let beams = [
        beam(10.0, 3e20, 0.75),
        beam(14.0, 2e21, 0.5),
        beam(26.0, 7e20, 1.0),
        beam(60.0, 1e21, 0.4),
    ];
    // In particular, the thin case must not compute absorption as 1-exp(-tau).
    for area in [1e14, 100.0, 1.0] {
        let result = evaluate_network_radiation_transport(
            &[NetworkRadiationCell::new(&state, &sources)],
            Area::new(area),
            &beams,
            NetworkRadiationSolve::new(controls(1e-6)),
            || false,
        )
        .unwrap();
        assert_eq!(result.cells()[0].initial_state(), &state);
        assert_eq!(result.cells()[0].sources(), sources);
        let mut expected_rates = [[0.0; 3]; 3];
        for (index, beam) in beams.iter().enumerate() {
            let energy_ev = beam.photons.photon_energy.value() / EV;
            let depths = channel_depths(&state, area, energy_ev);
            let depth = depths.iter().sum::<f64>();
            let mu = beam.direction_cosine.abs();
            let incoming = area * mu * beam.photons.number_flux_per_m2_s;
            let absorbed = incoming * -(-depth / mu).exp_m1();
            let outgoing = incoming * (-depth / mu).exp();
            for edge in 0..3 {
                let photons = if depth == 0.0 {
                    0.0
                } else {
                    absorbed * depths[edge] / depth
                };
                expected_rates[edge][0] += photons / NA;
                expected_rates[edge][1] += photons * energy_ev * EV;
                expected_rates[edge][2] += photons * (energy_ev - GAPS_EV[edge]).max(0.0) * EV;
            }
            let rate = result.beam_rates()[index];
            relative(rate.incident_photons_per_s, incoming, 2e-13);
            relative(rate.absorbed_photons_per_s, absorbed, 2e-12);
            relative(rate.transmitted_photons_per_s, outgoing, 2e-12);
            relative(rate.absorbed_power_w, absorbed * energy_ev * EV, 2e-12);
            assert!(rate.absorbed_photons_per_s <= incoming * (1.0 + 2e-15));
            assert!(rate.photon_residual_per_s.abs() <= 2e-15 * incoming);
            assert_eq!(result.cells()[0].beam_rates()[index], rate);
        }
        for (rate, expected) in result.cells()[0].channel_rates().iter().zip(expected_rates) {
            relative(rate.ionization_rate_mol_per_s, expected[0], 2e-12);
            relative(rate.absorbed_power_w, expected[1], 2e-12);
            relative(rate.electron_heating_power_w, expected[2], 2e-12);
        }
    }
}

#[test]
fn source_permutation_and_uniform_slab_subdivision_preserve_absorption() {
    let (state, sources) = fixture(false, [0.0; 3]);
    let beams = [beam(26.0, 2e21, 0.6), beam(60.0, 7e20, -0.3)];
    let solve = NetworkRadiationSolve::new(controls(1e-6));
    let whole = evaluate_network_radiation_transport(
        &[NetworkRadiationCell::new(&state, &sources)],
        Area::new(100.0),
        &beams,
        solve,
        || false,
    )
    .unwrap();
    let mut reversed = sources.clone();
    reversed.reverse();
    let reordered = evaluate_network_radiation_transport(
        &[NetworkRadiationCell::new(&state, &reversed)],
        Area::new(100.0),
        &beams,
        solve,
        || false,
    )
    .unwrap();
    assert_eq!(whole, reordered);

    let half = scaled(&state, 0.5);
    let split = evaluate_network_radiation_transport(
        &[
            NetworkRadiationCell::new(&half, &sources),
            NetworkRadiationCell::new(&half, &reversed),
        ],
        Area::new(100.0),
        &beams,
        solve,
        || false,
    )
    .unwrap();
    for (expected, actual) in whole.beam_rates().iter().zip(split.beam_rates()) {
        relative(
            actual.absorbed_photons_per_s,
            expected.absorbed_photons_per_s,
            2e-12,
        );
        relative(
            actual.transmitted_photons_per_s,
            expected.transmitted_photons_per_s,
            2e-12,
        );
    }
    for edge in 0..3 {
        let expected = whole.cells()[0].channel_rates()[edge];
        let left = split.cells()[0].channel_rates()[edge];
        let right = split.cells()[1].channel_rates()[edge];
        relative(
            left.ionization_rate_mol_per_s + right.ionization_rate_mol_per_s,
            expected.ionization_rate_mol_per_s,
            2e-12,
        );
        relative(
            left.absorbed_power_w + right.absorbed_power_w,
            expected.absorbed_power_w,
            2e-12,
        );
        relative(
            left.electron_heating_power_w + right.electron_heating_power_w,
            expected.electron_heating_power_w,
            2e-12,
        );
    }
}

#[test]
fn oblique_beams_enter_opposite_faces_and_shield_downstream_cells() {
    let (state, sources) = fixture(false, [0.0; 3]);
    let states = [scaled(&state, 0.25), scaled(&state, 0.75)];
    let beams = [beam(26.0, 2e21, 0.3), beam(60.0, 8e20, -0.8)];
    let area = 50.0;
    let result = evaluate_network_radiation_transport(
        &states
            .iter()
            .map(|state| NetworkRadiationCell::new(state, &sources))
            .collect::<Vec<_>>(),
        Area::new(area),
        &beams,
        NetworkRadiationSolve::new(controls(1e-6)),
        || false,
    )
    .unwrap();
    for (index, beam) in beams.iter().enumerate() {
        let energy_ev = beam.photons.photon_energy.value() / EV;
        let mu = beam.direction_cosine.abs();
        let boundary = area * mu * beam.photons.number_flux_per_m2_s;
        let mut incoming = boundary;
        let order = if beam.direction_cosine > 0.0 {
            [0, 1]
        } else {
            [1, 0]
        };
        for cell in order {
            let depth = channel_depths(&states[cell], area, energy_ev)
                .iter()
                .sum::<f64>()
                / mu;
            let absorbed = incoming * -(-depth).exp_m1();
            let outgoing = incoming * (-depth).exp();
            let actual = result.cells()[cell].beam_rates()[index];
            relative(actual.incident_photons_per_s, incoming, 2e-12);
            relative(actual.absorbed_photons_per_s, absorbed, 2e-12);
            relative(actual.transmitted_photons_per_s, outgoing, 2e-12);
            incoming = outgoing;
        }
        let total = result.beam_rates()[index];
        relative(total.incident_photons_per_s, boundary, 2e-12);
        relative(total.absorbed_photons_per_s, boundary - incoming, 2e-12);
        relative(total.transmitted_photons_per_s, incoming, 2e-12);
        assert!(result.cells()[order[1]].beam_rates()[index].incident_photons_per_s < boundary);
    }
}
