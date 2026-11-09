//! Original-source local-temperature ternary interfaces and atomic source gates.
//! Reuses the original pure-liquid Cp sources and mixture fixtures from the parent.

use super::*;
use fs_flux::miscibility::multicomponent::phase_field::{
    MulticomponentPhaseField, MulticomponentPhaseFieldNetwork,
    nonisothermal::{
        MulticomponentNonisothermalPhaseField,
        adaptive::MulticomponentNonisothermalPhaseFieldTimeControl,
        material::{
            MulticomponentPhaseFieldMaterialNonisothermalError as PhaseThermalError,
            ResolvedMulticomponentNonisothermalPhaseField,
        },
    },
};

fn phase_model(inverse_temperature: [f64; 9]) -> MulticomponentRegularSolutionV1 {
    MulticomponentRegularSolutionV1::new(
        (0..3).map(species).collect(),
        (0..3)
            .flat_map(|i| (0..3).map(move |j| if i == j { 0.0 } else { 4.0 }))
            .collect(),
        inverse_temperature.into(),
        Temperature::new(290.0),
        Temperature::new(600.0),
    )
    .unwrap()
}

fn phase_field(
    model: MulticomponentRegularSolutionV1,
    initial: [[f64; 3]; 2],
    kappa_j_m: f64,
) -> MulticomponentPhaseField {
    let transport = MaxwellStefanV1::new(
        [2, 0, 1]
            .map(|i| (species(i), MolarMass::new(MASSES[i])))
            .into(),
    )
    .unwrap();
    MulticomponentPhaseFieldNetwork::new(
        model,
        transport,
        Temperature::new(300.0),
        AmountConcentration::new(1000.0),
        vec![MulticomponentMiscibilityCell { volume_m3: 0.001 }; 2],
        vec![MulticomponentMiscibilityFace {
            cells: [0, 1],
            area_m2: 1.0,
            distance_m: 1.0,
            diffusivity_matrix: 0,
        }],
        vec![vec![0.0, 1e-6, 1e-4, 1e-6, 0.0, 1e-4, 1e-4, 1e-4, 0.0]],
        kappa_j_m,
        16 * 1024 * 1024,
        || false,
    )
    .unwrap()
    .bind(initial.map(Vec::from).into(), || false)
    .unwrap()
}

fn phase_time() -> MulticomponentNonisothermalPhaseFieldTimeControl {
    MulticomponentNonisothermalPhaseFieldTimeControl {
        thermal: thermal_time(),
        absolute_interfacial_energy_tolerance_j: 1e-8,
        relative_interfacial_energy_tolerance: 1e-5,
    }
}

// Independent integral of authored pure Cp, original B/T excess enthalpy,
// and the explicit kappa*A/(4*d)*sum(delta_x^2) interface functional. No caloric total,
// interface-energy getter or conservation receipt enters this oracle.
fn total_energy(field: &MulticomponentNonisothermalPhaseField) -> f64 {
    let bulk: f64 = field
        .temperatures()
        .iter()
        .enumerate()
        .map(|(cell, t)| {
            let n = field.cell_component_moles(cell).unwrap();
            let dt = t.value() - 300.0;
            let total: f64 = n.iter().sum();
            let sensible: f64 = (0..3)
                .map(|i| n[i] * (CP[i] * dt + 0.5 * CP_SLOPE[i] * dt * dt))
                .sum();
            let excess: f64 = (0..3)
                .flat_map(|i| ((i + 1)..3).map(move |j| R * B[3 * i + j] * n[i] * n[j] / total))
                .sum();
            sensible + excess
        })
        .sum();
    let interface: f64 = field
        .network()
        .faces()
        .iter()
        .map(|face| {
            let a = field.cell_mole_fractions(face.cells[0]).unwrap();
            let b = field.cell_mole_fractions(face.cells[1]).unwrap();
            0.25 * field.network().gradient_coefficient_j_m()
                * face.area_m2
                / face.distance_m
                * a.iter().zip(b).map(|(a, b)| (a - b) * (a - b)).sum::<f64>()
        })
        .sum();
    bulk + interface
}

fn assert_original_endpoint(
    field: &ResolvedMulticomponentNonisothermalPhaseField<'_, '_, '_>,
    source: &MulticomponentLiquidSolutionLaw<'_>,
    card: &MaterialCard,
) {
    assert!(std::ptr::eq(field.source(), source));
    assert_eq!(field.maximum_resolved_cells(), 2);
    assert_eq!(field.maximum_relative_volume_mismatch(), 1e-12);
    for (cell, resolved) in field.cells().iter().enumerate() {
        let t = field.field().temperatures()[cell];
        let n = field.field().cell_component_moles(cell).unwrap();
        let x = field.field().cell_mole_fractions(cell).unwrap();
        let expected = source.homogeneous(t, n, || false).unwrap();
        assert_eq!(resolved, &expected);
        assert_eq!(resolved.component_moles(), n);
        assert_eq!(resolved.thermodynamics().temperature(), t);
        assert_eq!(resolved.thermodynamics().model(), source.thermodynamics());
        near(resolved.volume_m3(), 0.001);
        near(resolved.molar_concentration_mol_m3(), 1000.0);
        let material = resolved.liquid().material();
        assert_eq!(material.card_identity(), card.content_hash());
        assert_eq!(material.state_coordinate("T"), Some(t.value()));
        assert_eq!(material.state_coordinate("pressure"), Some(100_000.0));
        near(material.state_coordinate("x_a").unwrap(), x[0]);
        near(material.state_coordinate("x_b").unwrap(), x[1]);
        let transport = resolved.transport_state().unwrap();
        near(transport.density_kg_m3(), 44.0 - 42.0 * x[0] - 16.0 * x[1]);
        near(
            transport.specific_heat_j_kg_k(),
            3500.0 + 500.0 * (t.value() - 300.0) / 300.0,
        );
        near(
            transport.thermal_conductivity_w_m_k(),
            0.3 + 0.3 * (t.value() - 300.0) / 300.0,
        );
        near(
            transport.kinematic_viscosity_m2_s(),
            1e-6 + 2e-6 * (t.value() - 300.0) / 300.0,
        );
        for property in material.properties() {
            let original = source
                .reference_state()
                .material()
                .property(property.requirement().name())
                .unwrap();
            assert_eq!(
                property.answer().receipt.selected,
                original.answer().receipt.selected
            );
        }
    }
}

#[test]
fn thermal_interfaces_retain_exact_calorics_and_refresh_actual_local_sources() {
    let card = temperature_source_card();
    let model = phase_model(B);
    let liquid = source(&card, model.clone(), MASSES);
    let caloric = caloric(&liquid);
    let numerical = phase_field(
        model.clone(),
        [[0.3, 0.3, 0.4], [0.31, 0.1, 0.59]],
        100.0,
    );
    let original = numerical
        .resolve_liquid_field(&liquid, 2, 1e-12, || false)
        .unwrap();
    let local = original
        .bind_nonisothermal(
            &caloric,
            vec![Temperature::new(320.0), Temperature::new(500.0)],
            faces(),
            bracket(),
            16 * 1024 * 1024,
            || false,
        )
        .unwrap();
    assert_original_endpoint(&local, &liquid, &card);
    assert!(std::ptr::eq(local.caloric_law(), &caloric));
    assert_eq!(local.field().caloric_model(), caloric.model());
    assert_eq!(local.network(), numerical.network());
    assert_eq!(local.field().transport_faces(), faces());
    // Matching values from an independently created liquid do not grant identity.
    let other_liquid = source(&card, model, MASSES);
    let other_caloric = self::caloric(&other_liquid);
    assert!(matches!(
        original.bind_nonisothermal(
            &other_caloric,
            local.field().temperatures().to_vec(),
            faces(),
            bracket(),
            16 * 1024 * 1024,
            || false,
        ),
        Err(PhaseThermalError::Invalid(_))
    ));
    let before = local.field().clone();
    let before_cells = local.cells().to_vec();
    let energy = total_energy(&before);
    // Absolute source integrals incur rounding; the physical allowance remains
    // the same closed-interval energy budget used by the numerical owner.
    let oracle_roundoff_j = 64.0 * f64::EPSILON * energy.abs().max(1.0);
    let controls = MulticomponentNonisothermalSolve::default();
    let raw = before.step_candidate(0.001, controls, || false).unwrap();
    let mut completed_polls = 0;
    let candidate = local
        .step_candidate(0.001, controls, || {
            completed_polls += 1;
            false
        })
        .unwrap();
    let end = candidate.field();
    assert_eq!(end.field(), raw.field());
    assert_eq!(candidate.step(), raw.step());
    assert_original_endpoint(end, &liquid, &card);
    assert!(std::ptr::eq(end.caloric_law(), &caloric));
    assert_eq!(end.field().elapsed_s(), 0.001);
    assert_eq!(end.network(), numerical.network());
    assert_eq!(end.field().transport_faces(), faces());
    assert_ne!(end.field().temperatures(), before.temperatures());
    assert!(
        (total_energy(end.field()) - energy).abs()
            <= controls.energy_tolerance_j + oracle_roundoff_j
    );
    assert!(candidate.step().energy_residual_j.abs() <= controls.energy_tolerance_j);
    assert!(candidate.step().entropy_change_j_per_k >= -controls.entropy_tolerance_j_per_k);
    assert!(
        (candidate.step().homogeneous_enthalpy_change_j
            + candidate.step().interfacial_energy_change_j)
            .abs()
            <= controls.energy_tolerance_j
    );
    for species in 0..3 {
        near(
            end.field().cell_component_moles(0).unwrap()[species]
                + end.field().cell_component_moles(1).unwrap()[species],
            before.cell_component_moles(0).unwrap()[species]
                + before.cell_component_moles(1).unwrap()[species],
        );
    }
    for stop in [1, completed_polls] {
        let mut polls = 0;
        assert!(matches!(
            local.step_candidate(0.001, controls, || {
                polls += 1;
                polls == stop
            }),
            Err(PhaseThermalError::Numerical(
                MulticomponentNonisothermalError::Cancelled
            ))
        ));
        assert_eq!(local.field(), &before);
        assert_eq!(local.cells(), before_cells.as_slice());
    }
    let replay = local.step_candidate(0.001, controls, || false).unwrap();
    assert_eq!(replay.field().field(), end.field());
    assert_eq!(replay.field().cells(), end.cells());
    assert_eq!(replay.step(), candidate.step());
    let (accepted, step) = candidate.into_parts();
    assert_eq!(step.duration_s, 0.001);
    assert!(std::ptr::eq(accepted.caloric_law(), &caloric));
    let raw = accepted
        .field()
        .advance_candidate(0.01, phase_time(), controls, || false)
        .unwrap();
    let advanced = accepted
        .advance_candidate(0.01, phase_time(), controls, || false)
        .unwrap();
    assert_eq!(advanced.field().field(), raw.field());
    assert_eq!(advanced.advance(), raw.advance());
    assert_original_endpoint(advanced.field(), &liquid, &card);
    // One accepted fixed step followed by one complete adaptive interval.
    assert!(
        (total_energy(advanced.field().field()) - energy).abs()
            <= 2.0 * controls.energy_tolerance_j + oracle_roundoff_j
    );
    let (advanced, report) = advanced.into_parts();
    assert!(report.accepted_intervals > 0);
    assert!(report.maximum_accepted_error_ratio <= 1.0);
    assert!(std::ptr::eq(advanced.caloric_law(), &caloric));
    assert!(std::ptr::eq(advanced.source(), &liquid));
    assert_eq!(local.field(), &before);
    assert_eq!(local.cells(), before_cells.as_slice());
}

#[test]
fn completed_thermal_interfaces_cannot_escape_original_source_or_volume_support() {
    let controls = MulticomponentNonisothermalSolve::default();
    for (initial, upper, perturbation) in [
        ([[0.31, 0.29, 0.4], [0.29, 0.31, 0.4]], 0.31, 0.0),
        ([[0.3, 0.3, 0.4], [0.31, 0.1, 0.59]], 1.0, 1.0),
    ] {
        let card = source_card(upper, perturbation);
        // A purely entropic nonconvex bulk law isolates phase enrichment from
        // excess mixing enthalpy. The energetic interface still exchanges heat.
        let model = phase_model([0.0; 9]);
        let liquid = source(&card, model.clone(), MASSES);
        let caloric = caloric(&liquid);
        let numerical = phase_field(model, initial, 0.01);
        let original = numerical
            .resolve_liquid_field(&liquid, 2, 1e-12, || false)
            .unwrap();
        let local = original
            .bind_nonisothermal(
                &caloric,
                vec![Temperature::new(350.0); 2],
                vec![MulticomponentNonisothermalFace {
                    thermal_conductivity_w_m_k: 0.0,
                    heat_of_transport_j_per_mol: vec![0.0; 3],
                }],
                bracket(),
                16 * 1024 * 1024,
                || false,
            )
            .unwrap();
        let before = local.field().clone();
        let before_cells = local.cells().to_vec();
        let raw = before.step_candidate(0.1, controls, || false).unwrap();
        if perturbation == 0.0 {
            assert!(raw.field().cell_mole_fractions(0).unwrap()[0] > upper);
        }
        before
            .advance_candidate(0.1, phase_time(), controls, || false)
            .unwrap();
        for error in [
            local.step_candidate(0.1, controls, || false).unwrap_err(),
            local
                .advance_candidate(0.1, phase_time(), controls, || false)
                .unwrap_err(),
        ] {
            if perturbation == 0.0 {
                assert!(matches!(
                    error,
                    PhaseThermalError::Material(Error::Source { cell: 0, .. })
                ));
            } else {
                assert!(matches!(
                    error,
                    PhaseThermalError::Material(Error::VolumeMismatch { .. })
                ));
            }
        }
        assert_eq!(local.field(), &before);
        assert_eq!(local.cells(), before_cells.as_slice());
        assert!(std::ptr::eq(local.source(), &liquid));
        assert!(std::ptr::eq(local.caloric_law(), &caloric));
    }
}
