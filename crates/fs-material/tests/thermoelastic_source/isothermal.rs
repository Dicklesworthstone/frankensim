//! G1/G3: original coefficient derivatives, reference inertia and strain domain.
use super::*;
use fs_blake3::ContentHash;
use fs_material::elastic::thermal::IsothermalThermoelasticInput;

fn inventory() -> IsothermalThermoelasticInput {
    IsothermalThermoelasticInput {
        reference_density_kg_m3: 6.0,
        mechanical_material_identity: ContentHash([61; 32]),
    }
}

fn matvec(c: &[[f64; 6]; 6], x: [f64; 6]) -> [f64; 6] {
    std::array::from_fn(|i| c[i].iter().zip(x).map(|(c, x)| c * x).sum())
}
fn dot(a: [f64; 6], b: [f64; 6]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}

#[test]
fn exact_isothermal_tensor_and_reference_inventory_survive_temperature_changes() {
    let card = source_card(-2.0, 0.0001, 0.0008, 0.000002, false);
    let law = law(&card);
    let reference = inventory();
    let mut previous_identity = None;
    for t in [281.0, 300.0, 321.0, 359.0] {
        let at = law.at_temperature(t).unwrap();
        let frozen = at.isothermal_elasticity(reference).unwrap();
        assert_eq!(frozen.source(), &at);
        assert_eq!(*frozen.input(), reference);
        assert_eq!(frozen.acoustic().density_kg_m3(), reference.reference_density_kg_m3);
        assert!(frozen.acoustic().is_explicitly_isotropic());
        for (actual, expected) in frozen.acoustic().stiffness_mandel_pa().iter().flatten()
            .zip(at.stiffness_mandel_pa().into_iter().flatten())
        {
            assert_eq!(actual.to_bits(), expected.to_bits());
        }
        let e = 1000.0 - 2.0 * (t - 300.0);
        let nu = 0.2 + 0.0008 * (t - 300.0);
        let (bulk, shear) = frozen.acoustic().isotropic_bulk_shear_pa().unwrap();
        close(bulk, e / (3.0 * (1.0 - 2.0 * nu)), 1e-12);
        close(shear, e / (2.0 * (1.0 + nu)), 1e-12);
        close(at.free_strain().unwrap(), 0.0001 * (t - 300.0)
            + 0.000001 * (t - 300.0).powi(2), 2e-17);
        let identity = frozen.acoustic().material_state_identity();
        assert_eq!(identity, at.isothermal_elasticity(reference).unwrap()
            .acoustic().material_state_identity());
        if let Some(previous) = previous_identity { assert_ne!(identity, previous); }
        previous_identity = Some(identity);
    }
}

#[test]
fn coefficient_derivatives_match_temperature_differences_and_the_native_potential() {
    let total = [0.013, -0.004, 0.002, 0.007, -0.003, 0.005];
    for parameters in [
        [-2.0, -0.0007, 0.0001, -0.000002],
        [2.0, 0.0008, -0.0001, 0.000003],
        [0.0, 0.001, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.000002],
        [0.0; 4],
    ] {
        let [et, nut, alpha, alphat] = parameters;
        let card = source_card(et, alpha, nut, alphat, false);
        let law = law(&card);
        let t = 321.0;
        let h = 0.05;
        let samples = [-2.0, -1.0, 0.0, 1.0, 2.0]
            .map(|step| law.at_temperature(t + step * h).unwrap());
        let at = &samples[2];
        let derivatives = at.temperature_derivatives().unwrap();
        let first = |values: [f64; 5]| {
            (values[0] - 8.0 * values[1] + 8.0 * values[3] - values[4]) / (12.0 * h)
        };
        let second = |values: [f64; 5]| {
            (-values[0] + 16.0 * values[1] - 30.0 * values[2]
                + 16.0 * values[3] - values[4]) / (12.0 * h * h)
        };
        for i in 0..6 {
            let free = std::array::from_fn(|k| samples[k].free_strain_mandel()[i]);
            close(derivatives.free_strain_mandel_per_k[i], first(free), 2e-13);
            close(derivatives.free_strain_mandel_per_k2[i], second(free), 2e-11);
            for j in 0..6 {
                let c = std::array::from_fn(|k| samples[k].stiffness_mandel_pa()[i][j]);
                close(derivatives.stiffness_mandel_pa_k[i][j], first(c), 2e-8);
                close(derivatives.stiffness_mandel_pa_k2[i][j], second(c), 2e-6);
            }
        }
        let moduli = samples.each_ref().map(|sample| sample.isothermal_elasticity(inventory())
            .unwrap().acoustic().isotropic_bulk_shear_pa().unwrap());
        let (kt, gt) = derivatives.isotropic_moduli_pa_k.unwrap();
        let (ktt, gtt) = derivatives.isotropic_moduli_pa_k2.unwrap();
        close(kt, first(moduli.map(|(k, _)| k)), 2e-8);
        close(gt, first(moduli.map(|(_, g)| g)), 2e-8);
        close(ktt, second(moduli.map(|(k, _)| k)), 2e-6);
        close(gtt, second(moduli.map(|(_, g)| g)), 2e-6);

        // Independent tensor work contractions check the public coefficient
        // derivatives against the original entropy/capacity/mixed response.
        let volume = 0.7;
        let potential = at.state(total, volume).unwrap();
        let e = potential.elastic_strain_mandel;
        let sigma = matvec(&at.stiffness_mandel_pa(), e);
        let first_stress = matvec(&derivatives.stiffness_mandel_pa_k, e);
        let second_stress = matvec(&derivatives.stiffness_mandel_pa_k2, e);
        let expansion = derivatives.free_strain_mandel_per_k;
        let expansion_stress = matvec(&at.stiffness_mandel_pa(), expansion);
        for i in 0..6 {
            close(potential.stress_temperature_derivative_mandel_pa_k[i],
                first_stress[i] - expansion_stress[i], 2e-11);
        }
        let entropy = volume * (dot(expansion, sigma) - 0.5 * dot(e, first_stress));
        let capacity = -t * volume * (0.5 * dot(e, second_stress)
            - 2.0 * dot(expansion, first_stress) + dot(expansion, expansion_stress)
            - dot(derivatives.free_strain_mandel_per_k2, sigma));
        close(potential.entropy_j_k, entropy, 2e-13);
        close(potential.heat_capacity_j_k, capacity, 2e-12);
    }
}

#[test]
fn frozen_view_preserves_the_total_and_elastic_component_domains() {
    let card = source_card(0.0, 0.0001, 0.0, 0.0, false);
    let law = law(&card);
    let cold = law.at_temperature(300.0).unwrap();
    let frozen = cold.isothermal_elasticity(inventory()).unwrap();
    assert_eq!(cold.maximum_absolute_mandel_strain(), 0.05);
    assert_eq!(frozen.acoustic().maximum_absolute_mandel_strain(), 0.05);
    assert_eq!(frozen.acoustic().frobenius_strain_limit(), None);
    // Every component is admitted although the Frobenius norm exceeds 0.05.
    let inside = [0.04; 6];
    assert!(frozen.source().state(inside, 1.0).is_ok());
    assert!(frozen.acoustic().admits_mandel_strain(&inside));
    assert!(!frozen.acoustic().admits_mandel_strain(&[0.051; 6]));
    let hot = law.at_temperature(360.0).unwrap().isothermal_elasticity(inventory()).unwrap();
    let total = [0.053, 0.0, 0.0, 0.0, 0.0, 0.0];
    let free = hot.source().free_strain_mandel();
    let elastic = std::array::from_fn(|i| total[i] - free[i]);
    assert!(hot.acoustic().admits_mandel_strain(&elastic));
    assert!(hot.source().state(total, 1.0).is_err());
    let total = [-0.048, 0.0, 0.0, 0.0, 0.0, 0.0];
    assert!(hot.source().state(total, 1.0).is_err());
}

#[test]
fn inventory_identity_and_invalid_reference_inputs_remain_explicit() {
    let card = source_card(-2.0, 0.0001, 0.0008, 0.000002, false);
    let at = law(&card).at_temperature(321.0).unwrap();
    let original = at.isothermal_elasticity(inventory()).unwrap();
    for input in [
        IsothermalThermoelasticInput { reference_density_kg_m3: 7.0, ..inventory() },
        IsothermalThermoelasticInput { mechanical_material_identity: ContentHash([62; 32]), ..inventory() },
    ] {
        let other = at.isothermal_elasticity(input).unwrap();
        assert_ne!(original.acoustic().material_state_identity(), other.acoustic().material_state_identity());
        assert_eq!(original.source(), other.source());
        assert_eq!(original.acoustic().stiffness_mandel_pa(), other.acoustic().stiffness_mandel_pa());
    }
    for density in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(at.isothermal_elasticity(IsothermalThermoelasticInput {
            reference_density_kg_m3: density, ..inventory()
        }).is_err());
    }
    assert!(at.isothermal_elasticity(IsothermalThermoelasticInput {
        mechanical_material_identity: ContentHash([0; 32]), ..inventory()
    }).is_err());
}
