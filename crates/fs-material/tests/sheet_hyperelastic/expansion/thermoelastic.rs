//! Original source paths, relaxed thickness and Green strain share one energy.
use super::*;

fn prepared(card: &MaterialCard, t: f64, isotropic: bool, exact: bool) -> ResolvedSheet {
    sheet(card, exact)
        .with_hyperelastic_free_expansion(
            card,
            path(card, "T", isotropic),
            t,
            if isotropic {
                SheetExpansionFrame::Isotropic
            } else {
                frame([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]])
            },
        )
        .unwrap()
}

#[test]
fn g1_source_membrane_entropy_capacity_and_thermal_stress_match_relaxed_energy() {
    let card = thermal_card("T", [290.0, 310.0]);
    let strain = [0.18, -0.04, 0.12];
    let temperature = 304.0;
    for isotropic in [true, false] {
        let sheet = prepared(&card, temperature, isotropic, false);
        let law = bound(&sheet);
        let state = law.thermoelastic_state(strain).unwrap();
        let mechanical = law.response(strain).unwrap();
        close(
            state.response().energy_j_m2,
            mechanical.energy_j_m2,
            1e-9,
            1e-11,
        );
        for i in 0..3 {
            close(
                state.response().resultants_n_m[i],
                mechanical.resultants_n_m[i],
                1e-8,
                1e-10,
            );
        }
        for i in 0..9 {
            close(
                state.response().tangent_n_m[i],
                mechanical.tangent_n_m[i],
                2e-7,
                2e-10,
            );
        }
        // Independent endpoint preparation and the ORIGINAL mechanical sheet
        // response; no thermal derivative enters this finite-difference oracle.
        let step = 0.01;
        let samples = [-2.0, -1.0, 1.0, 2.0].map(|offset| {
            let endpoint = sheet
                .at_hyperelastic_temperature(&card, temperature + offset * step)
                .unwrap();
            bound(&endpoint).response(strain).unwrap()
        });
        let derivative = |values: [f64; 4]| {
            (values[0] - 8.0 * values[1] + 8.0 * values[2] - values[3]) / (12.0 * step)
        };
        let energy_t = derivative(samples.map(|sample| sample.energy_j_m2));
        close(state.entropy_j_m2_k(), -energy_t, 2e-8, 3e-8);
        let energy_tt = (-samples[3].energy_j_m2 + 16.0 * samples[2].energy_j_m2
            - 30.0 * mechanical.energy_j_m2
            + 16.0 * samples[1].energy_j_m2
            - samples[0].energy_j_m2)
            / (12.0 * step * step);
        close(
            state.mechanical_heat_capacity_j_m2_k(),
            -temperature * energy_tt,
            2e-5,
            2e-5,
        );
        close(
            state.internal_energy_j_m2(),
            mechanical.energy_j_m2 + temperature * state.entropy_j_m2_k(),
            1e-12,
            1e-13,
        );
        for i in 0..3 {
            close(
                state.resultant_temperature_derivative_n_m_k()[i],
                derivative(samples.map(|sample| sample.resultants_n_m[i])),
                2e-8,
                2e-8,
            );
        }
        let fixed_full_f_capacity = law.reference_thickness_m()
            * state
                .native()
                .full_state()
                .mechanical_heat_capacity_j_m3_k();
        assert!(state.mechanical_heat_capacity_j_m2_k() > fixed_full_f_capacity);
        assert!((state.mechanical_heat_capacity_j_m2_k() - fixed_full_f_capacity).abs() > 1e-3);
        assert_eq!(state.source().material(), sheet.material());
        assert_eq!(state.strain(), &strain);
    }
}

#[test]
fn g1_entropy_constrained_green_tangent_uses_total_area_capacity_once() {
    let card = thermal_card("T", [290.0, 310.0]);
    let t = 304.0;
    let sheet = prepared(&card, t, false, false);
    let law = bound(&sheet);
    let strain = [0.18, 0.07, 0.09];
    let state = law.thermoelastic_state(strain).unwrap();
    // Authored caloric reference for this algebraic fixture, in J/(m² K).
    let capacity = 336.0;
    let adiabatic = state.adiabatic_response(capacity).unwrap();
    close(
        adiabatic.fixed_strain_heat_capacity_j_m2_k,
        capacity + state.mechanical_heat_capacity_j_m2_k(),
        0.0,
        1e-14,
    );
    let step = 1e-6;
    for column in 0..3 {
        let responses = [-1.0, 1.0].map(|sign| {
            let mut e = strain;
            e[column] += sign * step;
            let temperature = t + sign * step * adiabatic.temperature_gradient_k[column];
            let endpoint = sheet
                .at_hyperelastic_temperature(&card, temperature)
                .unwrap();
            let endpoint_law = bound(&endpoint);
            let endpoint_state = endpoint_law.thermoelastic_state(e).unwrap();
            let entropy = capacity * (temperature / t).ln() + endpoint_state.entropy_j_m2_k();
            (endpoint_state.response().resultants_n_m, entropy)
        });
        close(
            (responses[1].1 - responses[0].1) / (2.0 * step),
            0.0,
            2e-6,
            0.0,
        );
        for row in 0..3 {
            close(
                adiabatic.tangent_n_m[3 * row + column],
                (responses[1].0[row] - responses[0].0[row]) / (2.0 * step),
                2e-5,
                2e-7,
            );
        }
    }
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(state.adiabatic_response(bad).is_err());
    }
}

#[test]
fn g0_source_smoothness_and_constrained_thermal_scope_are_explicit() {
    let card = thermal_card("T", [290.0, 310.0]);
    let path = path_at(&card, "T", thermal_point("T", 299.0), false);
    // Preparation inserts the original datum as a knot. Its unchanged source
    // slope remains differentiable; the genuine measured kink at 300 K does not.
    let regular = path.at_temperature_with_smooth_tangent(299.0).unwrap();
    assert_eq!(regular.geometry(), path.at_temperature(299.0).unwrap());
    for bad in [290.0, 300.0, 310.0, f64::NAN, f64::INFINITY] {
        assert!(path.at_temperature_with_smooth_tangent(bad).is_err());
    }
    assert!(
        bound(&prepared(&card, 304.0, true, true))
            .thermoelastic_state([0.1, 0.0, 0.0])
            .is_err()
    );
    assert!(
        bound(&sheet(&card, false))
            .thermoelastic_state([0.1, 0.0, 0.0])
            .is_err()
    );
    assert!(
        bound(&prepared(&card, 300.0, true, false))
            .thermoelastic_state([0.1, 0.0, 0.0])
            .is_err()
    );
    let limited = thermal_card("T", [299.0, 304.0]);
    assert!(
        bound(&prepared(&limited, 304.0, true, false))
            .thermoelastic_state([0.1, 0.0, 0.0])
            .is_err()
    );
    let law = bound(&prepared(&card, 304.0, false, false));
    let jet = endpoint(&law).thermal_reference_jet().unwrap();
    let slopes = endpoint(&law)
        .source()
        .at_temperature_with_smooth_tangent(304.0)
        .unwrap()
        .principal_stretch_derivatives_per_k();
    assert_eq!(
        jet.first_derivative_per_k(),
        &[
            slopes[1], 0.0, 0.0, 0.0, slopes[0], 0.0, 0.0, 0.0, slopes[2]
        ]
    );
    assert_eq!(jet.second_derivative_per_k2(), &[0.0; 9]);
}

#[test]
fn g3_source_thermodynamics_retains_nonproportional_viscous_memory() {
    let card = thermal_card("T", [290.0, 310.0]);
    let t = 304.0;
    let sheet = prepared(&card, t, false, false);
    let source = bound(&sheet);
    let law = ViscoHyperelastic::new(
        source.model().law().clone(),
        vec![FiniteStrainBranch {
            energy: BranchEnergy::NeoHookean { mu_pa: 1.6e5 },
            tau_s: 0.2,
        }],
    )
    .unwrap()
    .with_reference_deformation(*endpoint(&source).reference_deformation());
    let strain = [0.07, 0.16, -0.06];
    let mut history = ViscoMembraneState::initial(&law, VolumetricTreatment::OriginalMaterial);
    for e in [[0.18, -0.02, 0.12], strain] {
        history = FrozenViscoMembrane::new(&law, &history, 0.08, sheet.thickness_m(), controls())
            .unwrap()
            .step_response(e)
            .unwrap()
            .state;
    }
    let before = history.clone();
    let thermal = source
        .frozen_visco_thermoelastic_state(&law, &history, strain)
        .unwrap();
    let native = FrozenViscoMembrane::new(&law, &history, 0.0, sheet.thickness_m(), controls())
        .unwrap()
        .response(strain)
        .unwrap();
    close(
        thermal.response().energy_j_m2,
        native.energy_j_m2,
        1e-8,
        1e-10,
    );
    for i in 0..3 {
        close(
            thermal.response().resultants_n_m[i],
            native.resultants_n_m[i],
            2e-7,
            2e-9,
        );
    }
    let ViscoMembraneState::Full(accepted) = &history else {
        unreachable!()
    };
    assert!(core::ptr::eq(
        thermal.viscoelastic_history().unwrap().1,
        accepted
    ));
    let step = 0.005;
    let energies = [-1.0, 1.0].map(|sign| {
        let target_sheet = sheet
            .at_hyperelastic_temperature(&card, t + sign * step)
            .unwrap();
        let target_source = bound(&target_sheet);
        let target = law
            .clone()
            .with_reference_deformation(*endpoint(&target_source).reference_deformation());
        // Native finite reference binding changes G while preserving Cv.
        // Its signed ledger is not used as an entropy or derivative oracle.
        let rebound =
            ViscoMembraneState::Full(target.change_reference_from(&law, accepted).unwrap().state);
        FrozenViscoMembrane::new(&target, &rebound, 0.0, sheet.thickness_m(), controls())
            .unwrap()
            .response(strain)
            .unwrap()
            .energy_j_m2
    });
    close(
        thermal.entropy_j_m2_k(),
        -(energies[1] - energies[0]) / (2.0 * step),
        2e-7,
        2e-6,
    );
    let virgin = ViscoMembraneState::initial(&law, VolumetricTreatment::OriginalMaterial);
    let erased = source
        .frozen_visco_thermoelastic_state(&law, &virgin, strain)
        .unwrap();
    assert!((erased.entropy_j_m2_k() - thermal.entropy_j_m2_k()).abs() > 1e-3);
    let foreign = law.clone().with_reference_deformation(
        fs_material::hyper::reference::ReferenceDeformation::identity(),
    );
    assert!(
        source
            .frozen_visco_thermoelastic_state(&foreign, &history, strain)
            .is_err()
    );
    assert_eq!(history, before);
}
