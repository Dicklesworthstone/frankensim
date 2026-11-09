//! One coupled gas trajectory against an independent two-cell ODE reference.
use super::*;
use fs_evidence::ValidityDomain;
use fs_lbm::conjugate::{AdaptiveSpeciesEnthalpyDiffusionControls, GasSpeciesDiffusivity};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::diffusion::{
    BINARY_DIFFUSIVITY_DIMS, BINARY_DIFFUSIVITY_PROPERTY, BinaryDiffusionMaterial,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

#[test]
fn g1_adaptive_nasa_enthalpy_and_binary_diffusivity_follow_the_coupled_ode() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let mass = [1.0, 2.0];
    let accepted = carrier(&domain, &mass);
    let initial = [[0.75, 0.25], [0.25, 1.75]];
    let a7 = -4000.0;
    let bank = source::moist_model(0.5, a7);
    let h = [
        enthalpy(&bank, initial[0], 320.0),
        enthalpy(&bank, initial[1], 480.0),
    ];
    let names = ids();
    let context = source::point(400.0, 100_000.0);
    let gate = CancelGate::new_clock_free();
    let typed = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(250.0, 0.1), (650.0, 0.26)],
                dims: BINARY_DIFFUSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    typed(QuantityKind::AbsoluteTemperature),
                    250.0,
                    650.0,
                )
                .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1e7),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "analytic binary D(T)=0.0004*T coupling fixture; no empirical gas claim"
                    .into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let material = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured/water binary pair".into(),
            phase: "gas".into(),
            process: "analytic adaptive transport fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let binding = BinaryDiffusionMaterial {
        species: names.clone(),
        material: &material,
        context: &context,
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    let calorics = GasSpeciesEnthalpySource::new_with_diffusivity(
        &bank,
        &names,
        &context,
        GasSpeciesDiffusivity::Binary(&binding),
        energy(),
        &gate,
    )
    .unwrap();
    let duration = 3.0;

    // Independent continuum ODE: q is A transferred left -> right; e is
    // energy received by the left cell. B counterdiffuses by exactly -q.
    // Integrate actual physical fluxes, not the production exponential sweeps
    // or NASA inversion. The fixture's Cp polynomial has an analytic inverse.
    let inventories = |q: f64| {
        [
            [initial[0][0] - q, initial[0][1] + q],
            [initial[1][0] + q, initial[1][1] - q],
        ]
    };
    let rhs = |[q, e]: [f64; 2]| {
        let species = inventories(q);
        let tl = temperature(mass[0] * h[0] + e, species[0], a7);
        let tr = temperature(mass[1] * h[1] - e, species[1], a7);
        // Unit cell volumes/area and half-cell distances of 0.5 m.
        // The native carrier densities are 1 and 2 kg/m³, not EOS replacements.
        let conductance = 1.0 / (0.5 / (mass[0] * 0.0004 * tl) + 0.5 / (mass[1] * 0.0004 * tr));
        let rate = conductance * (species[0][0] / mass[0] - species[1][0] / mass[1]);
        let face_h = partial_h(0.5 * tl + 0.5 * tr, a7);
        [rate, -rate * (face_h[0] - face_h[1])]
    };
    let reference = |steps: usize| {
        let dt = duration / steps as f64;
        let mut state = [0.0; 2];
        for _ in 0..steps {
            let k1 = rhs(state);
            let k2 = rhs(core::array::from_fn(|k| state[k] + 0.5 * dt * k1[k]));
            let k3 = rhs(core::array::from_fn(|k| state[k] + 0.5 * dt * k2[k]));
            let k4 = rhs(core::array::from_fn(|k| state[k] + dt * k3[k]));
            for k in 0..2 {
                state[k] += dt * (k1[k] + 2.0 * k2[k] + 2.0 * k3[k] + k4[k]) / 6.0;
            }
        }
        state
    };
    let oracle = reference(4096);
    let coarse_oracle = reference(2048);
    near(oracle[0], coarse_oracle[0], 1e-10);
    near(oracle[1], coarse_oracle[1], 1e-5);

    let limits = AdaptiveSpeciesEnthalpyDiffusionControls {
        solver: controls(),
        initial_substeps: 1,
        maximum_refinements: 10,
        maximum_source_evaluations: 1 << 18,
        absolute_temperature_tolerance_k: 2e-6,
        absolute_specific_enthalpy_tolerance_j_kg: 0.002,
        absolute_mass_fraction_tolerance: [1e-9; 2],
        relative_mass_fraction_tolerance: 1e-8,
        maximum_enthalpy_defect_j: 1e-4,
    };
    let before = accepted.clone();
    let result = accepted
        .diffuse_species_with_enthalpy_adaptive(
            &domain, &initial, &h, duration, &calorics, limits, &gate,
        )
        .unwrap();
    let expected_species = inventories(oracle[0]);
    for (c, sign) in [(0, 1.0), (1, -1.0)] {
        for k in 0..2 {
            near(result.species_mass_kg[c][k], expected_species[c][k], 1e-8);
        }
        let expected_h = h[c] + sign * oracle[1] / mass[c];
        let expected_t = temperature(mass[c] * expected_h, expected_species[c], a7);
        near(result.specific_enthalpy_j_kg[c], expected_h, 0.005);
        near(result.temperature_k[c], expected_t, 5e-6);
        near(
            result.diffusive_enthalpy_change_j[c],
            sign * oracle[1],
            0.005,
        );
    }
    // Refreshing D matters here: its initial-left value is 0.128 m²/s,
    // while the coupled trajectory raises it by more than 0.003 m²/s.
    // The accepted original source image is retained, with no extra inversion.
    let endpoint = &result.state;
    assert_eq!(endpoint.temperatures_k(), result.temperature_k.as_slice());
    assert_eq!(endpoint.source().source_packs(), bank.source_packs());
    assert!(endpoint.diffusivity_m2_s()[0] - 0.128 > 0.003);
    for c in 0..2 {
        near(
            endpoint.diffusivity_m2_s()[c],
            0.0004 * result.temperature_k[c],
            1e-12,
        );
    }
    assert!(result.report.refinements > 1);
    assert!(result.report.maximum_error_ratio <= 1.0);
    assert!(result.report.absolute_enthalpy_defect_j <= limits.maximum_enthalpy_defect_j);
    near(result.diffusive_enthalpy_change_j.iter().sum(), 0.0, 1e-6);
    assert_eq!(accepted, before);
}
