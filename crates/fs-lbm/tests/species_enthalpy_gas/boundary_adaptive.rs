//! Source-owned open gas diffusion against an independent scalar flux ODE.
use super::*;
use fs_evidence::ValidityDomain;
use fs_lbm::conjugate::{
    AdaptiveCoupledBoundarySpeciesEnthalpyDiffusionControls,
    AdaptiveSpeciesEnthalpyDiffusionControls, GasSpeciesDiffusivity,
    GasSpeciesEnthalpyBoundarySource, SpeciesDiffusionBoundary,
};
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
fn adaptive_original_gas_reservoir_refreshes_diffusivity_and_closes_external_energy() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let accepted = carrier(&domain, &[1.0]);
    let masses = [[0.9, 0.1]];
    let a7 = -4000.0;
    let bank = source::moist_model(0.5, a7);
    let h = [enthalpy(&bank, masses[0], 300.0)];
    let names = ids();
    let context = source::point(300.0, 100_000.0);
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
                source: "analytic binary D(T)=0.0004*T reservoir fixture; no empirical gas claim"
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
            process: "analytic open gas transport fixture".into(),
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
    let boundaries = [SpeciesDiffusionBoundary {
        cell: 0,
        face: 0,
        reservoir_mass_fractions: [0.1, 0.9],
        film_resistance_m2_s_kg: 2.5,
    }];
    let face_t = [600.0];
    let bound =
        GasSpeciesEnthalpyBoundarySource::new(&calorics, &boundaries, &face_t, &gate).unwrap();
    let duration = 3.0;
    let face_h = partial_h(face_t[0], a7);
    let enthalpy_difference = face_h[0] - face_h[1];

    // q is signed A received by the cell; B receives -q. One prescribed
    // local face T gives the exact energy invariant Q=q*(h_A-h_B).
    // Integrate the continuum flux directly, independently of the production
    // exponential midpoint maps, material interpolation and NASA inversion.
    let state_temperature =
        |q: f64| temperature(h[0] + q * enthalpy_difference, [0.9 + q, 0.1 - q], a7);
    let rhs = |q: f64| {
        // Unit carrier density/area/volume, half-cell width 0.5 m.
        let d = 0.0004 * state_temperature(q);
        let conductance = 1.0 / (0.5 / d + 2.5);
        conductance * (0.1 - (0.9 + q))
    };
    let reference = |steps: usize| {
        let dt = duration / steps as f64;
        let mut q = 0.0;
        for _ in 0..steps {
            let k1 = rhs(q);
            let k2 = rhs(q + 0.5 * dt * k1);
            let k3 = rhs(q + 0.5 * dt * k2);
            let k4 = rhs(q + dt * k3);
            q += dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        }
        q
    };
    let oracle = reference(4096);
    near(oracle, reference(2048), 1e-12);
    let expected_heat = oracle * enthalpy_difference;
    let expected_t = state_temperature(oracle);
    let limits = || AdaptiveCoupledBoundarySpeciesEnthalpyDiffusionControls {
        integrator: AdaptiveSpeciesEnthalpyDiffusionControls {
            solver: controls(),
            initial_substeps: 1,
            maximum_refinements: 10,
            maximum_source_evaluations: 1 << 18,
            absolute_temperature_tolerance_k: 2e-6,
            absolute_specific_enthalpy_tolerance_j_kg: 0.002,
            absolute_mass_fraction_tolerance: [1e-9; 2],
            relative_mass_fraction_tolerance: 1e-8,
            maximum_enthalpy_defect_j: 1e-4,
        },
        absolute_boundary_transfer_tolerance_kg: [1e-8; 2],
        absolute_boundary_enthalpy_tolerance_j: 0.002,
    };
    let before = (accepted.clone(), masses, h);
    let result = accepted
        .diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive(
            &domain,
            &masses,
            &h,
            &boundaries,
            duration,
            &bound,
            limits(),
            &gate,
        )
        .unwrap();
    for (k, dm) in [(0, oracle), (1, -oracle)] {
        near(result.step.species_mass_kg[0][k], masses[0][k] + dm, 1e-8);
        near(result.boundary.mass_into_domain_kg[k], dm, 1e-8);
        near(result.boundary.faces[0].mass_into_cell_kg[k], dm, 1e-8);
        assert!(
            result.estimated_boundary_transfer_error_kg[k]
                <= limits().absolute_boundary_transfer_tolerance_kg[k]
        );
    }
    near(
        result.step.specific_enthalpy_j_kg[0],
        h[0] + expected_heat,
        0.005,
    );
    near(
        result.step.diffusive_enthalpy_change_j[0],
        expected_heat,
        0.005,
    );
    near(
        result.thermal_boundary.enthalpy_into_domain_j,
        expected_heat,
        0.005,
    );
    near(
        result.thermal_boundary.faces[0].enthalpy_into_cell_j,
        expected_heat,
        0.005,
    );
    near(result.step.temperature_k[0], expected_t, 5e-6);
    let endpoint = &result.step.state;
    assert_eq!(
        endpoint.temperatures_k(),
        result.step.temperature_k.as_slice()
    );
    assert_eq!(endpoint.species(), &names);
    assert_eq!(endpoint.source().source_packs(), bank.source_packs());
    near(
        endpoint.diffusivity_m2_s()[0],
        0.0004 * result.step.temperature_k[0],
        1e-12,
    );
    near(
        endpoint.transport(0).unwrap().specific_heat_j_kg_k(),
        capacity(result.step.species_mass_kg[0]) * (3.0 + 0.001 * expected_t),
        1e-5,
    );
    assert!(endpoint.diffusivity_m2_s()[0] - 0.12 > 0.016);
    // The unrefreshed D=0.12 trajectory misses more than 0.01 kg of A.
    let frozen_conductance: f64 = 1.0 / (0.5 / 0.12 + 2.5);
    let frozen_q = -0.8 * (-(-frozen_conductance * duration).exp_m1());
    assert!((oracle - frozen_q).abs() > 0.01);
    assert!(result.step.report.refinements > 1);
    assert!(result.step.report.maximum_error_ratio <= 1.0);
    assert!(
        result.step.report.absolute_enthalpy_defect_j
            <= limits().integrator.maximum_enthalpy_defect_j
    );
    assert!(
        result.estimated_boundary_enthalpy_error_j
            <= limits().absolute_boundary_enthalpy_tolerance_j
    );

    // An otherwise valid native boundary cannot borrow a trace bound to a
    // different film/composition/face descriptor through the generic driver.
    let mismatched = [SpeciesDiffusionBoundary {
        film_resistance_m2_s_kg: 3.0,
        ..boundaries[0]
    }];
    assert!(matches!(
        accepted.diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive(
            &domain, &masses, &h, &mismatched, duration, &bound, limits(), &gate,
        ),
        Err(ChtError::InvalidInput { reason, .. })
            if reason.contains("bound face-temperature descriptor")
    ));
    assert_eq!(accepted, before.0);
    assert_eq!(masses, before.1);
    assert_eq!(h, before.2);
}
