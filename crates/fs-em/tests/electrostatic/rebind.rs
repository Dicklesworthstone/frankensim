//! G1/G3/G4: native material substitution at retained charge, with independent
//! capacitor/work and relaxation oracles. Coefficients are synthetic SI data.
use super::*;
use fs_em::{
    conduction::PlanarConduction,
    electrostatic::relaxation::{
        DielectricMaterialRebindControls, DielectricRelaxationControls, DielectricRelaxationError,
        PlanarLeakyDielectric,
    },
};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    field_storage::planar::{PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES, PlanarDielectricLaw},
};

fn conductor(
    epsilon: [f64; 2],
    sigma: f64,
    basis: ElectricalPropertyBasis,
    revision: u32,
) -> ResolvedConductorState {
    let electrical_value = match basis {
        ElectricalPropertyBasis::Conductivity => sigma,
        ElectricalPropertyBasis::Resistivity => sigma.recip(),
    };
    let properties = [
        (
            ELECTRIC_PERMITTIVITY_PROPERTY,
            ELECTRIC_PERMITTIVITY_DIMS,
            epsilon[0],
        ),
        (
            PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES[0],
            ELECTRIC_PERMITTIVITY_DIMS,
            epsilon[0],
        ),
        (
            PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES[1],
            ELECTRIC_PERMITTIVITY_DIMS,
            epsilon[1],
        ),
        (basis.property_name(), basis.dimensions(), electrical_value),
    ];
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic material-replacement dielectric".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("authored dielectric {revision}"),
            phase: "solid".into(),
            process: "analytical fixture".into(),
            revision,
        },
        claims,
        vec![],
    )
    .unwrap();
    let material = resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ResolvedConductorState::try_from_material_with_basis(material, basis).unwrap()
}
fn evolution() -> DielectricRelaxationControls {
    DielectricRelaxationControls {
        linear: config(),
        charge_tolerance_c: 1e-10,
        current_tolerance_a: 1e-10,
        energy_tolerance_j: 1e-9,
    }
}
fn replacement() -> DielectricMaterialRebindControls {
    DielectricMaterialRebindControls {
        linear: config(),
        charge_tolerance_c: 1e-10,
        energy_tolerance_j: 1e-9,
        maximum_absolute_work_j: 100.,
    }
}
fn field(held: &[usize], offset: f64) -> PlanarLeakyDielectric {
    let source = conductor([2.; 2], 0.5, ElectricalPropertyBasis::Conductivity, 0);
    let conduction = PlanarConduction::new(
        layered(0.5, 1., 0.2),
        vec![source; 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    PlanarLeakyDielectric::new(
        conduction,
        [0., 1.5, 3., 0., 1.5, 3.].map(|v| v + offset).to_vec(),
        held,
        7.,
        limits(),
        evolution(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn free_charge_survives_material_and_basis_swap_then_drives_new_relaxation() {
    let old = field(&[0, 3], 0.);
    let source = conductor([4.; 2], 2., ElectricalPropertyBasis::Resistivity, 1);
    let (new, work) = old
        .rebind_material_candidate(
            vec![source.clone(); 4],
            vec![PlanarDielectricLaw::Isotropic; 4],
            limits(),
            replacement(),
            &mut || false,
        )
        .unwrap();
    // A=0.2, L=1: C doubles from 0.4 to 0.8 F. At conserved Q=1.2 C,
    // V halves and U=Q^2/(2C) halves. Ground exchanges no work.
    near(new.voltage_v()[2], 1.5);
    near(new.voltage_v()[1], 0.75);
    near(old.stored_energy_j(), 1.8);
    near(new.stored_energy_j(), 0.9);
    near(work.parameter_work_j, -0.9);
    assert!(work.supplied_work_j.abs() < 1e-12);
    assert!(work.energy_residual_j.abs() < 1e-12);
    assert_eq!(new.time_s(), old.time_s());
    for v in [1, 2, 4, 5] {
        assert_eq!(
            new.nodal_charge_c()[v].to_bits(),
            old.nodal_charge_c()[v].to_bits()
        );
    }
    for i in 0..4 {
        assert_ne!(
            work.material_state_identities[i][0],
            work.material_state_identities[i][1]
        );
        assert_eq!(new.conduction().conductors()[i], source);
    }
    let next = new
        .step_to_candidate(
            7.5,
            vec![source; 4],
            &[(0, 0.), (3, 0.)],
            limits(),
            evolution(),
            &mut || false,
        )
        .unwrap();
    // New tau=epsilon/sigma=2 s, midpoint decay over dt=0.5 is 7/9.
    near(next.field().voltage_v()[2], 1.5 * 7. / 9.);
    near(next.step().heat_j, 0.9 * (1. - (7.0_f64 / 9.).powi(2)));
    assert!(next.step().parameter_work_j.abs() < 1e-12);
    assert!(next.step().energy_residual_j.abs() < 1e-12);
    // A true new card remains distinct from ordinary pinned-source refresh.
    assert!(
        old.step_to_candidate(
            7.5,
            new.conduction().conductors().to_vec(),
            &[(0, 0.), (3, 0.)],
            limits(),
            evolution(),
            &mut || false,
        )
        .is_err()
    );
}

#[test]
fn held_electrode_work_and_material_work_close_separately_and_are_gauge_invariant() {
    for offset in [0., 1024.] {
        let old = field(&[0, 3, 2, 5], offset);
        let (new, work) = old
            .rebind_material_candidate(
                vec![conductor([4.; 2], 2., ElectricalPropertyBasis::Conductivity, 1); 4],
                vec![PlanarDielectricLaw::Isotropic; 4],
                limits(),
                replacement(),
                &mut || false,
            )
            .unwrap();
        // Fixed V=3: delta_U=delta_C*V^2/2=1.8 J, source work=V*delta_Q=3.6 J.
        // The replacement owner receives 1.8 J; it is not irreversible heat.
        near(work.supplied_work_j, 3.6);
        near(work.parameter_work_j, -1.8);
        near(new.stored_energy_j(), 3.6);
        near(
            work.nodal_electrode_charge_c[2].unwrap() + work.nodal_electrode_charge_c[5].unwrap(),
            1.2,
        );
        assert!(work.electrode_charge_residual_c.abs() < 1e-10);
        assert!(work.energy_residual_j.abs() < 1e-9);
        for v in [0, 3, 2, 5] {
            assert_eq!(new.voltage_v()[v], old.voltage_v()[v]);
        }
    }
}

#[test]
fn heterogeneous_replacement_redistributes_field_without_erasing_interface_charge() {
    let old = field(&[0, 3], 0.);
    let changed = conductor([4.; 2], 1., ElectricalPropertyBasis::Conductivity, 1);
    let same = old.conduction().conductors()[2].clone();
    let (new, work) = old
        .rebind_material_candidate(
            vec![changed.clone(), changed, same.clone(), same],
            vec![PlanarDielectricLaw::Isotropic; 4],
            limits(),
            replacement(),
            &mut || false,
        )
        .unwrap();
    // Series layers: C=A/(L_left/epsilon_left+L_right/epsilon_right)=8/15 F.
    // Original Q=1.2 C remains, so the left/right voltage drops are 0.75/1.5 V.
    near(new.voltage_v()[1], 0.75);
    near(new.voltage_v()[2], 2.25);
    near(new.stored_energy_j(), 1.35);
    near(work.parameter_work_j, -0.45);
    for v in [1, 2, 4, 5] {
        assert_eq!(
            new.nodal_charge_c()[v].to_bits(),
            old.nodal_charge_c()[v].to_bits()
        );
    }
    for d in new.displacement_c_m2(&mut || false).unwrap() {
        near(d[0], -6.);
        assert!(d[1].abs() < 1e-10);
    }
    // Undoing only the material data restores the physical field and returns
    // the opposite work; it never restores by discarding a previous state.
    let (restored, reverse) = new
        .rebind_material_candidate(
            old.conduction().conductors().to_vec(),
            vec![PlanarDielectricLaw::Isotropic; 4],
            limits(),
            replacement(),
            &mut || false,
        )
        .unwrap();
    near(restored.voltage_v()[2], 3.);
    near(reverse.parameter_work_j, 0.45);
}

#[test]
fn conductivity_only_replacement_preserves_history_exactly_and_changes_native_heat() {
    let old = field(&[0, 3], 0.);
    let mut controls = replacement();
    controls.maximum_absolute_work_j = 0.;
    let (new, work) = old
        .rebind_material_candidate(
            vec![conductor([2.; 2], 2., ElectricalPropertyBasis::Resistivity, 2); 4],
            vec![PlanarDielectricLaw::Isotropic; 4],
            limits(),
            controls,
            &mut || false,
        )
        .unwrap();
    assert_eq!(new.voltage_v(), old.voltage_v());
    assert_eq!(new.nodal_charge_c(), old.nodal_charge_c());
    assert_eq!(work.parameter_work_j, 0.);
    assert_eq!(work.supplied_work_j, 0.);
    assert!(work.linear.is_none());
    let a = old
        .step_to_candidate(
            7.5,
            old.conduction().conductors().to_vec(),
            &[(0, 0.), (3, 0.)],
            limits(),
            evolution(),
            &mut || false,
        )
        .unwrap();
    let b = new
        .step_to_candidate(
            7.5,
            new.conduction().conductors().to_vec(),
            &[(0, 0.), (3, 0.)],
            limits(),
            evolution(),
            &mut || false,
        )
        .unwrap();
    // Same field, new tau=1 s, compared with original tau=4 s.
    near(b.field().voltage_v()[2], 3. * 0.6);
    near(b.step().heat_j, 1.8 * (1. - 0.6_f64.powi(2)));
    assert!(b.step().heat_j > a.step().heat_j);
}

#[test]
fn rotated_dielectric_substitution_uses_tensor_cross_terms() {
    let old = field(&[0, 1, 2, 3, 4, 5], 0.);
    let law = super::anisotropic::law(core::f64::consts::FRAC_PI_4);
    let (new, work) = old
        .rebind_material_candidate(
            vec![conductor([2., 6.], 1., ElectricalPropertyBasis::Conductivity, 3); 4],
            vec![law; 4],
            limits(),
            replacement(),
            &mut || false,
        )
        .unwrap();
    // E=(-3,0), rotated epsilon=[[4,-2],[-2,4]], D=(-12,6).
    // Fixed potential: U=3.6 J, work_electrodes=3.6 J, work_material=-1.8 J.
    near(new.stored_energy_j(), 3.6);
    near(work.supplied_work_j, 3.6);
    near(work.parameter_work_j, -1.8);
    for d in new.displacement_c_m2(&mut || false).unwrap() {
        near(d[0], -12.);
        near(d[1], 6.);
    }
    assert!(new.dielectric().permittivity_f_m().is_none());
}

#[test]
fn work_refusal_and_late_cancellation_leave_original_state_available_for_retry() {
    let old = field(&[0, 3], 0.);
    let source = conductor([4.; 2], 2., ElectricalPropertyBasis::Resistivity, 4);
    let mut controls = replacement();
    controls.maximum_absolute_work_j = 0.;
    assert!(matches!(
        old.rebind_material_candidate(
            vec![source.clone(); 4],
            vec![PlanarDielectricLaw::Isotropic; 4],
            limits(),
            controls,
            &mut || false,
        ),
        Err(DielectricRelaxationError::Invalid(_))
    ));
    let mut polls = 0;
    let (_, expected) = old
        .rebind_material_candidate(
            vec![source.clone(); 4],
            vec![PlanarDielectricLaw::Isotropic; 4],
            limits(),
            replacement(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        old.rebind_material_candidate(
            vec![source.clone(); 4],
            vec![PlanarDielectricLaw::Isotropic; 4],
            limits(),
            replacement(),
            &mut || {
                calls += 1;
                calls == polls
            },
        ),
        Err(DielectricRelaxationError::Cancelled)
    ));
    let (retry, actual) = old
        .rebind_material_candidate(
            vec![source; 4],
            vec![PlanarDielectricLaw::Isotropic; 4],
            limits(),
            replacement(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(actual.parameter_work_j, expected.parameter_work_j);
    assert_eq!(old.voltage_v()[2], 3.);
    near(old.stored_energy_j(), 1.8);
    near(retry.voltage_v()[2], 1.5);
}
