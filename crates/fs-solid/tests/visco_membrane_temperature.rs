//! G1/G3/G4: spatial source temperatures drive retained native natural strain.

#[path = "visco_membrane/fixture.rs"]
mod fixture;
use fixture::*;

use std::sync::Arc;

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, PropertyClaim, PropertyKey,
    PropertyValue, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::sheet::hyperelastic::{HyperelasticMembrane, expansion::SheetExpansionFrame};
use fs_material::sheet::{MembraneLaw, ResolvedSheet, SheetCapabilities, resolve_sheet};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::transport::expansion::{
    FreeExpansionGeometryQuery, FreeExpansionStrainSources, PreparedFreeExpansion,
};
use fs_material::visco_hyper::membrane::ViscoMembraneState;
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_solid::membrane::MembraneError;
use fs_solid::membrane::inertia::MembraneMassImage;
use fs_solid::membrane::visco::dynamics::{FiniteViscoMembraneMotion, FiniteViscoMotionControls};
use fs_solid::membrane::visco::{
    FiniteViscoMembrane, FiniteViscoMembraneError, OwnedFiniteViscoMembraneMotion,
    ViscoelasticLawField,
};

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            temperature,
        )
        .unwrap()
}

fn thermal_source(exact: bool) -> (MaterialCard, ResolvedSheet) {
    let original = source_card();
    let validity = ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        290.0,
        310.0,
    );
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        claim.validity = validity.clone();
        claims.insert_claim(claim).unwrap();
    }
    let key = PropertyKey::new("thermal_engineering_length_change", Dims::NONE);
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                // Nonzero datum and a slope change expose compounding or use
                // of 1 + strain instead of the original length ratio.
                knots: vec![(290.0, -0.01), (300.0, 0.02), (305.0, 0.04), (310.0, 0.09)],
                dims: Dims::NONE,
            },
            interpolation: InterpolationPolicy::LinearInside,
            validity: validity.clone(),
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: original.models()[0].provenance.clone(),
        })
        .unwrap();
    let mut models = original.models().to_vec();
    models[0].validity = validity;
    let card = MaterialCard::assemble(original.id().clone(), claims, models).unwrap();
    let reference = point(300.0);
    let path = PreparedFreeExpansion::try_new(
        &card,
        FreeExpansionGeometryQuery {
            reference: &reference,
            strains: FreeExpansionStrainSources::Isotropic {
                strain: &key,
                selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            },
            temperature_axis: "T",
            lower: &point(290.0),
            upper: &point(310.0),
        },
    )
    .unwrap();
    let sheet = resolve_sheet(
        &card,
        &reference,
        SheetCapabilities {
            membrane: false,
            bending: false,
        },
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let sheet = if exact {
        sheet.with_incompressible_hyperelastic_membrane(&card, 0, local())
    } else {
        sheet.with_hyperelastic_membrane(&card, 0, local())
    }
    .unwrap()
    .with_hyperelastic_free_expansion(&card, path, 300.0, SheetExpansionFrame::Isotropic)
    .unwrap();
    (card, sheet)
}

fn bound(sheet: &ResolvedSheet) -> Arc<HyperelasticMembrane> {
    let MembraneLaw::Hyperelastic(law) = sheet.membrane_law().unwrap() else {
        unreachable!()
    };
    law
}

fn accepted(exact: bool) -> (MaterialCard, OwnedFiniteViscoMembraneMotion) {
    let (card, sheet) = thermal_source(exact);
    let law = material(&sheet);
    let laws: Vec<_> = (0..4)
        .map(|index| {
            let mut law = law.clone();
            law.branches[0].tau_s = 0.1 + 0.15 * index as f64;
            law
        })
        .collect();
    let patch = patch(diamond(), sheet)
        .with_reference_areal_density(
            vec![
                [0.7, 1.0, 1.4],
                [1.1, 0.8, 0.9],
                [0.6, 1.3, 1.2],
                [0.9, 0.95, 1.05],
            ],
            || false,
        )
        .unwrap();
    let membrane =
        FiniteViscoMembrane::try_new(&patch, ViscoelasticLawField::PerElement(&laws), || false)
            .unwrap();
    let motion = FiniteViscoMembraneMotion::try_new(
        membrane,
        &[[0.0; 3]; 5],
        MembraneMassImage::Consistent,
        || false,
    )
    .unwrap();
    (card, OwnedFiniteViscoMembraneMotion::from_native(&motion))
}

fn metrics(state: &ViscoMembraneState) -> &[[f64; 6]] {
    match state {
        ViscoMembraneState::Full(state) => &state.inverse_viscous_cauchy_green,
        ViscoMembraneState::Incompressible(state) => &state.inverse_viscous_cauchy_green,
    }
}

fn engineering_strain(temperature: f64) -> f64 {
    if temperature <= 300.0 {
        -0.01 + 0.003 * (temperature - 290.0)
    } else if temperature <= 305.0 {
        0.02 + 0.004 * (temperature - 300.0)
    } else {
        0.04 + 0.01 * (temperature - 305.0)
    }
}

#[test]
fn regional_source_heating_creates_constraint_stress_and_native_relaxation_heat() {
    for exact in [false, true] {
        let (card, owned) = accepted(exact);
        let before = owned.as_native();
        let temperatures = [310.0, 304.0, 298.0, 306.0];
        let event = owned
            .with_temperatures(&[&card; 4], &temperatures, controls(), || false)
            .unwrap();
        let view = event.motion.as_native();
        assert_eq!(
            view.membrane().displacement_m().as_ptr(),
            before.membrane().displacement_m().as_ptr()
        );
        assert_eq!(view.velocity_mps().as_ptr(), before.velocity_mps().as_ptr());
        assert_eq!(view.membrane().time_s(), 0.0);
        assert_eq!(view.mass_image(), before.mass_image());
        assert_eq!(
            event.motion.patch().mesh().nodes,
            owned.patch().mesh().nodes
        );
        assert_eq!(
            event.motion.patch().mesh().elems,
            owned.patch().mesh().elems
        );
        assert_eq!(event.motion.patch().grain_rad(), owned.patch().grain_rad());
        assert_eq!(
            event.motion.patch().reference_corner_areal_density_kg_m2(),
            owned.patch().reference_corner_areal_density_kg_m2()
        );
        assert_eq!(
            event.motion.patch().nodal_mass_kg(),
            owned.patch().nodal_mass_kg()
        );
        assert_eq!(event.energy.dissipated_j, 0.0);
        assert_eq!(event.impulse_kg_m_s, [0.0; 3]);
        assert_eq!(event.kinetic_energy_change_j, 0.0);
        assert!(event.energy.stored_after_j > event.energy.stored_before_j);
        assert!(
            view.membrane()
                .triangle_resultants_n_m()
                .iter()
                .flatten()
                .any(|x| x.abs() > 1e-5)
        );
        let ViscoelasticLawField::PerElement(laws) = event.motion.native_laws() else {
            unreachable!()
        };
        let ViscoelasticLawField::PerElement(original) = owned.native_laws() else {
            unreachable!()
        };
        for (index, temperature) in temperatures.into_iter().enumerate() {
            let sheet = &event.motion.patch().sheets()[index];
            let law = bound(sheet);
            let reference = law.free_expansion_reference().unwrap();
            let g = (1.0 + engineering_strain(temperature)) / 1.02;
            for stretch in law.reference_stretches() {
                close(stretch, g, 3e-14);
            }
            assert_eq!(reference.temperature_k(), temperature);
            assert_eq!(reference.reference_context(), &point(300.0));
            assert_eq!(reference.current_context(), &point(temperature));
            assert_eq!(
                reference.source(),
                bound(&owned.patch().sheets()[index])
                    .free_expansion_reference()
                    .unwrap()
                    .source()
            );
            assert_eq!(sheet.material(), owned.patch().sheets()[index].material());
            assert_eq!(
                laws[index].reference_deformation(),
                reference.reference_deformation()
            );
            assert_eq!(
                laws[index],
                original[index]
                    .clone()
                    .with_reference_deformation(*reference.reference_deformation())
            );
            let state = &view.membrane().triangle_states()[index];
            assert_eq!(
                metrics(state),
                metrics(&before.membrane().triangle_states()[index])
            );
            assert_eq!(state.dissipated_pa(), 0.0);
            if exact {
                close(state.deformation_gradient()[8], g * g * g, 3e-12);
                assert!(view.membrane().triangle_pressure_pa()[index].is_some());
            }
        }
        let zero = [[0.0; 3]; 5];
        let step = event
            .motion
            .advance_with_pressure(
                0.05,
                &zero,
                0.0,
                &prescribed(&zero),
                settings(),
                FiniteViscoMotionControls {
                    material: controls(),
                    maximum_numerical_kinetic_loss_j: 1e-12,
                },
                || false,
            )
            .unwrap();
        assert_eq!(step.motion.as_native().membrane().displacement_m(), &zero);
        assert_eq!(step.motion.as_native().velocity_mps(), &zero);
        assert_eq!(step.motion.as_native().membrane().time_s(), 0.05);
        assert!(step.material_energy.dissipated_j > 0.0);
        assert!(step.element_dissipated_j.iter().all(|heat| *heat > 0.0));
        close(
            step.element_dissipated_j.iter().sum(),
            step.material_energy.dissipated_j,
            1e-13,
        );
        assert!(
            step.equilibrium
                .residual_n
                .iter()
                .flatten()
                .any(|r| r.abs() > 1e-5)
        );
        assert_eq!(
            step.motion.patch().nodal_mass_kg(),
            owned.patch().nodal_mass_kg()
        );
    }
}

#[test]
fn temperature_revisits_use_the_original_engineering_length_datum() {
    for exact in [false, true] {
        let (card, owned) = accepted(exact);
        let hot = owned
            .with_temperatures(&[&card; 4], &[310.0; 4], controls(), || false)
            .unwrap();
        let target = [304.0, 291.0, 306.0, 309.0];
        let revisited = hot
            .motion
            .with_temperatures(&[&card; 4], &target, controls(), || false)
            .unwrap();
        let direct = owned
            .with_temperatures(&[&card; 4], &target, controls(), || false)
            .unwrap();
        assert_eq!(
            revisited.motion.patch().sheets(),
            direct.motion.patch().sheets()
        );
        let a = revisited.motion.as_native();
        let b = direct.motion.as_native();
        close(
            a.membrane().stored_energy_j(),
            b.membrane().stored_energy_j(),
            1e-11,
        );
        for (x, y) in a
            .membrane()
            .triangle_first_piola_pa()
            .iter()
            .flatten()
            .zip(b.membrane().triangle_first_piola_pa().iter().flatten())
        {
            close(*x, *y, 1e-9);
        }
        let returned = revisited
            .motion
            .with_temperatures(&[&card; 4], &[300.0; 4], controls(), || false)
            .unwrap();
        assert_eq!(returned.motion.patch().sheets(), owned.patch().sheets());
        close(
            returned.motion.as_native().membrane().stored_energy_j(),
            0.0,
            1e-11,
        );
        assert_eq!(returned.energy.dissipated_j, 0.0);
        assert_eq!(returned.motion.as_native().membrane().time_s(), 0.0);
        for state in returned.motion.as_native().membrane().triangle_states() {
            assert_eq!(state.dissipated_pa(), 0.0);
        }
    }
}

#[test]
fn missing_sources_coverage_work_and_final_cancellation_leave_the_owner_unchanged() {
    let (card, owned) = accepted(true);
    let before = owned.as_native();
    let foreign = source_card();
    let temperatures = [309.0; 4];
    for (cards, values) in [
        (vec![&card; 3], temperatures.to_vec()),
        (vec![&card; 4], vec![309.0; 3]),
        (vec![&card; 4], vec![309.0, 311.0, 309.0, 309.0]),
        (vec![&card; 4], vec![309.0, f64::NAN, 309.0, 309.0]),
        (vec![&card, &card, &foreign, &card], temperatures.to_vec()),
    ] {
        assert!(
            owned
                .with_temperatures(&cards, &values, controls(), || false)
                .is_err()
        );
    }
    assert!(matches!(
        owned.with_temperatures(
            &[&card, &card, &foreign, &card],
            &temperatures,
            controls(),
            || false
        ),
        Err(FiniteViscoMembraneError::Native(MembraneError::Material {
            element: 2,
            ..
        }))
    ));
    let plain = source(true);
    let plain_law = material(&plain);
    let plain_patch = patch(triangle(), plain);
    let plain_membrane = FiniteViscoMembrane::try_new(
        &plain_patch,
        ViscoelasticLawField::Uniform(&plain_law),
        || false,
    )
    .unwrap();
    let plain_motion = FiniteViscoMembraneMotion::try_new(
        plain_membrane,
        &[[0.0; 3]; 3],
        MembraneMassImage::Lumped,
        || false,
    )
    .unwrap();
    assert!(
        OwnedFiniteViscoMembraneMotion::from_native(&plain_motion)
            .with_temperatures(&[&foreign], &[305.0], controls(), || false)
            .is_err()
    );
    let mut polls = 0;
    let accepted = owned
        .with_temperatures(&[&card; 4], &temperatures, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(accepted.energy.absolute_endpoint_work_defect_j > 0.0);
    let mut strict = controls();
    strict.maximum_endpoint_work_defect_j = 0.0;
    assert!(matches!(
        owned.with_temperatures(&[&card; 4], &temperatures, strict, || false),
        Err(FiniteViscoMembraneError::EnergyGate { .. })
    ));
    let mut seen = 0;
    assert!(matches!(
        owned.with_temperatures(&[&card; 4], &temperatures, controls(), || {
            seen += 1;
            seen == polls
        }),
        Err(FiniteViscoMembraneError::Native(MembraneError::Cancelled))
    ));
    assert_eq!(seen, polls);
    let after = owned.as_native();
    assert_eq!(
        after.membrane().triangle_states(),
        before.membrane().triangle_states()
    );
    assert_eq!(
        after.membrane().triangle_states().as_ptr(),
        before.membrane().triangle_states().as_ptr()
    );
    assert_eq!(
        after.membrane().displacement_m().as_ptr(),
        before.membrane().displacement_m().as_ptr()
    );
    assert_eq!(
        after.velocity_mps().as_ptr(),
        before.velocity_mps().as_ptr()
    );
    assert!(owned.patch().sheets().iter().all(|sheet| {
        bound(sheet)
            .free_expansion_reference()
            .unwrap()
            .temperature_k()
            == 300.0
    }));
}
