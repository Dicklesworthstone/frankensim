//! Independent nonisothermal film dynamics and outgoing radiative work.
use super::*;
use fs_render::thermal_surface::heating::{
    CoatingOpticsModel, CoatingSurfaceSource, StatefulCoatingLayer,
    StatefulTransmittingSurfaceIlluminationSource, SurfaceCoatingLaw, SurfaceIncidence,
};
use std::sync::Arc;

const MASSES: [f64; 2] = [0.4, 1.2];
const INITIAL_H: [f64; 2] = [0.0, 8.0];
const AUTHORED_W: [f64; 2] = [0.6, 0.2];
const COSINE: f64 = 0.9;

fn varying_film(index: usize, temperature: f64) -> CoherentAbsorbingCoating {
    let delta = temperature - 300.0;
    if index == 0 {
        film(1.35 + 0.03 * delta, 0.04 + 0.025 * delta, 75.0)
    } else {
        film(2.3 - 0.02 * delta, 0.2 + 0.015 * delta, 120.0)
    }
}

// State = two native specific enthalpies, two reflected spectral energies,
// two transmitted spectral energies. Evaluate only the original optical work
// owner: no constitutive-table resolution, heat packet or observer is reused.
fn rate(state: [f64; 6]) -> [f64; 6] {
    let films = std::array::from_fn::<_, 2, _>(|i| varying_film(i, 300.0 + state[i] / 4.0));
    let layers = films.each_ref().map(CoatingLayer::Absorbing);
    let stack = DielectricCoatingStack::try_mixed(&layers).unwrap();
    let mut result = [0.0; 6];
    for lane in 0..2 {
        let response = stack
            .layer_power_balance(WAVELENGTHS[lane], true, COSINE, 1.5, 1.0)
            .unwrap();
        for (i, a) in response.layer_absorbed().iter().enumerate() {
            result[i] +=
                POWER[lane] * (POLARIZATION[lane] * a[0] + (1.0 - POLARIZATION[lane]) * a[1]);
        }
        let q = 2.0 * POLARIZATION[lane] - 1.0;
        result[2 + lane] = POWER[lane]
            * (response.fresnel.reflection.0[0][0] + q * response.fresnel.reflection.0[0][1]);
        result[4 + lane] = POWER[lane]
            * (response.fresnel.transmission.0[0][0] + q * response.fresnel.transmission.0[0][1]);
    }
    for i in 0..2 {
        result[i] = (result[i] + AUTHORED_W[i]) / MASSES[i];
    }
    result
}

fn rk4() -> [f64; 6] {
    let mut state = [INITIAL_H[0], INITIAL_H[1], 0.0, 0.0, 0.0, 0.0];
    let dt = 1.0 / 256.0;
    for _ in 0..256 {
        let k1 = rate(state);
        let k2 = rate(std::array::from_fn(|i| state[i] + dt * 0.5 * k1[i]));
        let k3 = rate(std::array::from_fn(|i| state[i] + dt * 0.5 * k2[i]));
        let k4 = rate(std::array::from_fn(|i| state[i] + dt * k3[i]));
        state = std::array::from_fn(|i| {
            state[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
        });
    }
    state
}

fn state_of(thermal: &DepositionSimulation<'_>, optical: &SurfaceOpticalHeating<2>) -> [f64; 6] {
    [
        thermal.cells()[&CellId(11)].specific_h_j_kg,
        thermal.cells()[&CellId(12)].specific_h_j_kg,
        optical.energy().reflected_j[0],
        optical.energy().reflected_j[1],
        optical.energy().transmitted_j[0],
        optical.energy().transmitted_j[1],
    ]
}

#[test]
fn independent_reverse_film_states_and_outgoing_radiation_converge_to_native_optical_ode() {
    let cells = (0..2)
        .map(|i| {
            let mut cell = cell(11 + i as u64, MASSES[i], 1.0, AUTHORED_W[i]);
            cell.specific_h_j_kg = INITIAL_H[i];
            cell
        })
        .collect();
    let native = thermal(cells, &SENSIBLE, vec![]);
    let source = StatefulTransmittingSurfaceIlluminationSource {
        source_id: 1,
        films: (0..2)
            .map(|i| StatefulCoatingLayer {
                receiver: receiver(11 + i as u64),
                optics: CoatingOpticsModel::Law(Arc::new(
                    SurfaceCoatingLaw::try_temperature(
                        optical_source(),
                        &[
                            (300.0, varying_film(i, 300.0)),
                            (330.0, varying_film(i, 330.0)),
                        ],
                        2,
                    )
                    .unwrap(),
                )),
            })
            .collect(),
        incidence: SurfaceIncidence {
            wavelengths_nm: WAVELENGTHS,
            incident_eta: [1.5; 2],
            incidence_cosine: COSINE,
            s_power_fraction: POLARIZATION,
        },
        transmitted_eta: [1.0; 2],
        reverse: true,
        source_identity: SOURCE,
        power: power(),
    };
    let mut limit = budget();
    limit.max_source_evaluations = 4096;
    let input = SurfaceOpticalHeating::try_new_coating_stacks(
        &[CoatingSurfaceSource::StatefulTransmitting(&source)],
        &native,
        limit,
        4,
        tolerance(),
        || Ok(()),
    )
    .unwrap();
    let initial = [INITIAL_H[0], INITIAL_H[1], 0.0, 0.0, 0.0, 0.0];
    let initial_rate = rate(initial);
    let (first, optical_first, _) = native
        .advance_observed(&input, 0.125, config(0.125), || Ok(()))
        .unwrap();
    let first_state = state_of(&first, &optical_first);
    for i in 0..6 {
        close(first_state[i], initial[i] + 0.125 * initial_rate[i]);
    }
    let (second, optical_second, _) = first
        .advance_observed(&optical_first, 0.25, config(0.125), || Ok(()))
        .unwrap();
    let second_state = state_of(&second, &optical_second);
    let second_rate = rate(first_state);
    for i in 0..6 {
        close(second_state[i], first_state[i] + 0.125 * second_rate[i]);
    }
    // Both receiving films evolve independently, but the coherent stack solves
    // once when refreshed. The outgoing transmission changes with their states.
    assert_eq!(optical_second.response_evaluations(), 1);
    assert!((second_rate[4] - initial_rate[4]).abs() > 1e-3);
    assert!((second_rate[5] - initial_rate[5]).abs() > 1e-3);
    let reference = rk4();
    let mut errors = Vec::new();
    for dt in [0.0625, 0.03125, 0.015625] {
        let (hot, optical, audit) = native
            .advance_observed(&input, 1.0, config(dt), || Ok(()))
            .unwrap();
        let result = state_of(&hot, &optical);
        let h_error: f64 = (0..2).map(|i| (result[i] - reference[i]).abs()).sum();
        let outgoing_error: f64 = (2..6).map(|i| (result[i] - reference[i]).abs()).sum();
        errors.push([h_error, outgoing_error]);
        assert_eq!(optical.response_evaluations(), (1.0 / dt) as usize - 1);
        close(
            hot.energy_j() - native.energy_j(),
            optical.energy().thermal_input_j + 0.8,
        );
        close(
            audit.external_heat_j,
            optical.energy().thermal_input_j + 0.8,
        );
        close(audit.residual_j, 0.0);
        for lane in 0..2 {
            assert!(optical.energy().transmitted_j[lane] > 0.0);
            close(optical.energy().incident_j[lane], POWER[lane]);
            close(
                optical.energy().reflected_j[lane]
                    + optical.energy().transmitted_j[lane]
                    + optical.energy().absorbed_j[lane],
                POWER[lane],
            );
        }
    }
    for component in 0..2 {
        assert!(
            errors[1][component] < 0.7 * errors[0][component],
            "{errors:?}"
        );
        assert!(
            errors[2][component] < 0.7 * errors[1][component],
            "{errors:?}"
        );
    }
    assert!(errors[2][0] < 0.1 && errors[2][1] < 0.02, "{errors:?}");
    assert_eq!(native.time_s(), 0.0);
    assert_eq!(input.time_s(), 0.0);
    assert_eq!(input.energy().thermal_input_j, 0.0);
}
