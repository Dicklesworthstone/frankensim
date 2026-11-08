use super::*;
use crate::chemistry::diffuse_recombination::{
    DiffuseRecombinationCell, DiffuseRecombinationSolve, evaluate_diffuse_recombination,
};
use crate::chemistry::ionization_network::transport::diffuse::*;

fn diffuse_controls(relative_tolerance: f64) -> NetworkDiffuseSolve {
    NetworkDiffuseSolve::new(
        NetworkRadiationSolve::new(controls(1e-6)),
        DiffuseRecombinationSolve {
            relative_tolerance,
            ..DiffuseRecombinationSolve::default()
        },
    )
}

fn mixture_state(
    original: &FrozenIonizedGasStateV1,
    amounts: [f64; 6],
    temperature_k: f64,
) -> FrozenIonizedGasStateV1 {
    FrozenIonizedGasStateV1::new(
        original.thermodynamics(),
        &amounts.map(Amount::new),
        original.volume(),
        original.heavy_temperature(),
        Temperature::new(temperature_k),
        1e-12,
        || false,
    )
    .unwrap()
}

// Independent infinite-domain constant-sigma Milne moments, specialized only
// to the existing manufactured H/He fixture. The original finite sources extend
// hundreds of thermal energies beyond these moments' support here.
fn milne(edge: usize, state: &FrozenIonizedGasStateV1) -> [f64; 3] {
    let theta = KB * state.electron_temperature().value();
    let gap = GAPS_EV[edge] * EV;
    let statistical_ratio = f64::from(WEIGHTS[edge][0]) / f64::from(WEIGHTS[edge][1]);
    let factor = statistical_ratio * (2.0 / core::f64::consts::PI).sqrt()
        / (ME.powf(1.5) * 299_792_458.0_f64.powi(2));
    let pairs = (NA * state.amounts()[EDGES[edge][1]].value()) * (NA * state.amounts()[5].value())
        / state.volume().value();
    let photons = factor * SIGMA[edge] / theta.sqrt()
        * (gap * gap + 2.0 * gap * theta + 2.0 * theta * theta)
        * pairs;
    let kinetic = factor
        * SIGMA[edge]
        * theta.sqrt()
        * (gap * gap + 4.0 * gap * theta + 6.0 * theta * theta)
        * pairs;
    [photons, gap * photons + kinetic, kinetic]
}

#[test]
fn one_source_recovers_existing_diffuse_transport_and_stable_escape_limits() {
    let (state, sources) = fixture(false, [0.0; 3]);
    let solve = diffuse_controls(1e-7);
    let exact = milne(0, &state);
    // Independent angular integrals already used by the single-source tests:
    // beta(tau)=int_0^1 (1-exp(-tau/mu))/(tau/mu) dmu.
    for (tau, escaped_fraction) in [
        (1e-6, 0.999_992_630_852_386_5),
        (1.0, 0.390_308_032_802_239_9),
        (1e20, 5e-21),
    ] {
        let area = Area::new(SIGMA[0] * NA * state.amounts()[0].value() / tau);
        let result = evaluate_network_diffuse_recombination(
            &[NetworkRadiationCell::new(&state, &sources[..1])],
            area,
            &[],
            solve,
            || false,
        )
        .unwrap();
        let existing = evaluate_diffuse_recombination(
            &[DiffuseRecombinationCell {
                state: &state,
                source: &sources[0],
            }],
            area,
            solve.radiation,
            || false,
        )
        .unwrap();
        let balance = result.balance();
        let channel = result.cells()[0].channel_rates()[0];
        relative(balance.emitted_photons_per_s, exact[0], 3e-7);
        relative(balance.emitted_power_w, exact[1], 3e-7);
        relative(channel.captured_electron_power_w, exact[2], 3e-7);
        relative(
            balance.escaped_photons_per_s,
            exact[0] * escaped_fraction,
            3e-6,
        );
        relative(balance.escaped_power_w, exact[1] * escaped_fraction, 3e-6);
        relative(
            balance.absorbed_photons_per_s,
            exact[0] * (1.0 - escaped_fraction),
            3e-6,
        );
        relative(-channel.net_gas_power_w, exact[1] * escaped_fraction, 3e-6);
        relative(
            -channel.net_ionization_rate_mol_per_s * NA,
            exact[0] * escaped_fraction,
            3e-6,
        );
        relative(
            -channel.net_electron_power_w,
            exact[2] * escaped_fraction,
            3e-6,
        );
        relative(
            balance.escaped_power_w,
            existing.balance().escaped_power_w,
            3e-6,
        );
        relative(
            channel.photoelectron_power_w,
            existing.cells()[0].balance().photoelectron_power_w,
            3e-6,
        );
        assert_eq!(channel.cross_channel_absorbed_photons_per_s, 0.0);
        assert_eq!(channel.nonlocal_absorbed_photons_per_s, 0.0);
        assert_eq!(result.cells()[0].initial_state(), &state);
        assert_eq!(result.cells()[0].sources(), &sources[..1]);
        relative(
            result.interfaces()[0].left_power_w,
            result.interfaces()[1].right_power_w,
            3e-13,
        );
        assert_eq!(result.interfaces()[0].right_photons_per_s, 0.0);
        assert_eq!(result.interfaces()[1].left_photons_per_s, 0.0);
    }
}

#[test]
fn helium_recombination_locally_ionizes_hydrogen_with_the_original_gap_difference() {
    let (original, sources) = fixture(true, [0.0; 3]);
    let state = mixture_state(&original, [1.0, 0.0, 0.5, 0.2, 0.0, 0.2], 20000.0);
    let exact = milne(1, &state);
    // Only He I recombines. All its photons exceed both retained thresholds,
    // so the two absorbers' constant opacity fractions are exactly 5/7 and 2/7.
    // No emitted spectrum, rate evaluator or production quadrature is used by
    // this reference. The trapped case must retain escape below the rounding
    // of the much larger cross-species chemical and electron-energy exchange.
    for (tau, beta) in [(1.0, 0.390_308_032_802_239_9), (1e20, 5e-21)] {
        let area = Area::new(NA * (SIGMA[0] + 0.5 * SIGMA[1]) / tau);
        let result = evaluate_network_diffuse_recombination(
            &[NetworkRadiationCell::new(&state, &sources[..2])],
            area,
            &[],
            diffuse_controls(1e-7),
            || false,
        )
        .unwrap();
        let cell = &result.cells()[0];
        let hydrogen = cell.channel_rates()[0];
        let helium = cell.channel_rates()[1];
        let h_fraction = (1.0 - beta) * 5.0 / 7.0;
        let he_fraction = (1.0 - beta) * 2.0 / 7.0;
        let h_absorption = exact[0] * h_fraction;
        let he_absorption = exact[0] * he_fraction;
        let gap_difference = (GAPS_EV[1] - GAPS_EV[0]) * EV;
        let h_heating = h_fraction * (exact[2] + gap_difference * exact[0]);
        relative(
            hydrogen.cross_channel_absorbed_photons_per_s,
            h_absorption,
            3e-6,
        );
        relative(helium.self_absorbed_photons_per_s, he_absorption, 3e-6);
        assert_eq!(hydrogen.self_absorbed_photons_per_s, 0.0);
        assert_eq!(hydrogen.emitted_photons_per_s, 0.0);
        assert_eq!(helium.cross_channel_absorbed_photons_per_s, 0.0);
        relative(hydrogen.photoelectron_power_w, h_heating, 3e-6);
        relative(helium.photoelectron_power_w, he_fraction * exact[2], 3e-6);
        relative(helium.captured_electron_power_w, exact[2], 3e-7);
        relative(
            hydrogen.net_ionization_rate_mol_per_s,
            h_absorption / NA,
            3e-6,
        );
        relative(
            helium.net_ionization_rate_mol_per_s,
            (he_absorption - exact[0]) / NA,
            3e-6,
        );
        relative(
            hydrogen.chemical_storage_rate_w,
            GAPS_EV[0] * EV * h_absorption,
            3e-6,
        );
        relative(
            helium.chemical_storage_rate_w,
            GAPS_EV[1] * EV * (he_absorption - exact[0]),
            3e-6,
        );
        relative(hydrogen.net_gas_power_w, h_fraction * exact[1], 3e-6);
        relative(helium.net_gas_power_w, (he_fraction - 1.0) * exact[1], 3e-6);
        relative(result.balance().escaped_power_w, beta * exact[1], 3e-6);
        relative(
            result.balance().escaped_photons_per_s,
            beta * exact[0],
            3e-6,
        );
        let direct = cell.balance();
        relative(-direct.net_gas_power_w, beta * exact[1], 3e-6);
        relative(
            -direct.net_ionization_rate_mol_per_s * NA,
            beta * exact[0],
            3e-6,
        );
        relative(
            direct.net_electron_power_w,
            h_fraction * gap_difference * exact[0] - beta * exact[2],
            3e-6,
        );
        assert!(direct.net_electron_power_w > 0.0);
        assert!(direct.net_gas_power_w < 0.0);
    }
}

#[test]
fn diffuse_competition_is_canonical_and_uniform_subdivision_preserves_the_field() {
    let (state, sources) = fixture(false, [0.0; 3]);
    let beams = [beam(14.0, 1e22, 0.6), beam(60.0, 2e21, -0.4)];
    let area = Area::new(100.0);
    let solve = diffuse_controls(1e-7);
    let whole = evaluate_network_diffuse_recombination(
        &[NetworkRadiationCell::new(&state, &sources)],
        area,
        &beams,
        solve,
        || false,
    )
    .unwrap();
    let mut reversed = sources.clone();
    reversed.reverse();
    let reordered = evaluate_network_diffuse_recombination(
        &[NetworkRadiationCell::new(&state, &reversed)],
        area,
        &beams,
        solve,
        || false,
    )
    .unwrap();
    assert_eq!(whole, reordered);
    let half = scaled(&state, 0.5);
    let halves = evaluate_network_diffuse_recombination(
        &[
            NetworkRadiationCell::new(&half, &sources),
            NetworkRadiationCell::new(&half, &reversed),
        ],
        area,
        &beams,
        solve,
        || false,
    )
    .unwrap();
    relative(
        halves.balance().emitted_photons_per_s,
        whole.balance().emitted_photons_per_s,
        3e-6,
    );
    relative(
        halves.balance().escaped_photons_per_s,
        whole.balance().escaped_photons_per_s,
        3e-6,
    );
    relative(
        halves.balance().escaped_power_w,
        whole.balance().escaped_power_w,
        3e-6,
    );
    relative(
        halves.balance().absorbed_power_w,
        whole.balance().absorbed_power_w,
        3e-6,
    );
    let refined_gas = halves
        .cells()
        .iter()
        .map(|c| c.balance().net_gas_power_w)
        .sum::<f64>();
    relative(
        refined_gas,
        whole.cells()[0].balance().net_gas_power_w,
        3e-6,
    );
    assert!(halves.cells().iter().any(|c| {
        c.channel_rates()
            .iter()
            .any(|r| r.nonlocal_absorbed_photons_per_s > 0.0)
    }));

    // Prescribed beams keep their original finite budget alongside recycled
    // photons; adding diffuse emission cannot consume an incident photon twice.
    let incident = evaluate_network_radiation_transport(
        &[NetworkRadiationCell::new(&state, &sources)],
        area,
        &beams,
        solve.network,
        || false,
    )
    .unwrap();
    assert_eq!(whole.beam_rates(), incident.beam_rates());
    relative(
        whole.cells()[0].balance().net_gas_power_w,
        whole.balance().incident_power_w - whole.balance().escaped_power_w,
        3e-6,
    );
}
