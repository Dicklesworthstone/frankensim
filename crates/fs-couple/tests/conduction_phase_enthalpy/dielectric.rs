//! Native spatial relaxation and its enthalpy feedback share the existing fixture.
use super::*;
#[path = "dielectric/anisotropic.rs"]
mod anisotropic;
use fs_em::electrostatic::relaxation::{DielectricRelaxationControls, PlanarLeakyDielectric};
use fs_material::{
    conductor::ResolvedConductorState,
    field_storage::{ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY},
};

fn material(slope: f64, upper_fraction: f64, varying_epsilon: bool) -> MaterialCard {
    let base = card(slope, upper_fraction);
    let mut claims = base.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(ELECTRIC_PERMITTIVITY_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS),
            value: if varying_epsilon {
                PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300., 1.), (350., 3.)],
                    dims: ELECTRIC_PERMITTIVITY_DIMS,
                }
            } else {
                PropertyValue::Scalar {
                    value: 2.,
                    dims: ELECTRIC_PERMITTIVITY_DIMS,
                }
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                300.,
                350.,
            ),
            interpolation: if varying_epsilon {
                InterpolationPolicy::LinearInside
            } else {
                InterpolationPolicy::ConstantWithinValidity
            },
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "authored nondispersive dielectric fixture; not empirical material data"
                    .into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}
fn conductor(card: &MaterialCard, context: &QueryPoint, tensor: bool) -> ResolvedPlanarConductor {
    let basis = ElectricalPropertyBasis::Resistivity;
    let names = if tensor {
        PRINCIPAL_RESISTIVITY_PROPERTIES.to_vec()
    } else {
        vec![basis.property_name()]
    };
    let mut requirements: Vec<_> = names
        .iter()
        .map(|&name| {
            ScalarPropertyRequirement::try_new(
                name,
                basis.dimensions(),
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()
        })
        .collect();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            ELECTRIC_PERMITTIVITY_PROPERTY,
            ELECTRIC_PERMITTIVITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let bundle = resolve_material_state_point(
        card,
        context,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    if tensor {
        let frame = PlanarElectricalFrame::new(
            ContentHash([1; 32]),
            ContentHash([2; 32]),
            core::f64::consts::FRAC_PI_4,
        )
        .unwrap();
        ResolvedOrthotropicConductor2::try_from_phase_material(bundle, basis, frame, FRACTION)
            .unwrap()
            .into()
    } else {
        ResolvedConductorState::try_from_material_with_basis(bundle, basis)
            .unwrap()
            .with_liquid_mass_fraction_axis(FRACTION)
            .unwrap()
            .into()
    }
}
fn electric_controls() -> DielectricRelaxationControls {
    DielectricRelaxationControls {
        linear: controls().linear,
        charge_tolerance_c: 1e-9,
        current_tolerance_a: 1e-9,
        energy_tolerance_j: 1e-8,
    }
}
fn field(
    card: &MaterialCard,
    context: &QueryPoint,
    tensor: bool,
    voltage: f64,
) -> PlanarLeakyDielectric<ResolvedPlanarConductor> {
    let source = conductor(card, context, tensor);
    let conduction =
        PlanarConduction::new(mesh(), vec![source; 4], limits(), &mut || false).unwrap();
    PlanarLeakyDielectric::new(
        conduction,
        vec![0., voltage, voltage, 0., voltage / 2.],
        &[0, 3],
        0.,
        limits(),
        electric_controls(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn native_dielectric_relaxation_matches_epsilon_over_sigma_and_local_physical_heat() {
    let material = material(0., 1., false);
    let context = point(310., 0.5);
    for tensor in [false, true] {
        let old = field(&material, &context, tensor, 3.);
        let next = old
            .step_to_candidate(
                0.5,
                old.conduction().conductors().to_vec(),
                &[(0, 0.), (3, 0.)],
                limits(),
                electric_controls(),
                &mut || false,
            )
            .unwrap();
        // Homogeneous charge relaxation has tau=epsilon/sigma=2 s. Midpoint
        // gives a=(2*tau-dt)/(2*tau+dt)=7/9; the actual field is affine.
        let a = 7. / 9.;
        let mid = 3. * (1. + a) / 2.;
        near(next.field().voltage_v()[1], 3. * a);
        near(next.field().voltage_v()[4], 1.5 * a);
        near(old.stored_energy_j(), 9.);
        near(next.field().stored_energy_j(), 9. * a * a);
        near(next.step().heat_j, 9. * (1. - a * a));
        near(next.step().supplied_work_j, 0.);
        for i in 0..4 {
            near(next.step().cell_power_w[i], mid * mid / 4.);
            for axis in 0..2 {
                near(
                    next.step().conduction_current_density_a_m2[i][axis],
                    -mid * core::f64::consts::FRAC_1_SQRT_2,
                );
                near(
                    next.step().conduction_current_density_a_m2[i][axis]
                        + next.step().displacement_current_density_a_m2[i][axis],
                    0.,
                );
            }
        }
        for (&before, &after) in old
            .nodal_charge_c()
            .iter()
            .zip(next.field().nodal_charge_c())
        {
            near(after, a * before);
        }
        let continuing = next
            .field()
            .step_to_candidate(
                1.,
                next.field().conduction().conductors().to_vec(),
                &[(0, 0.), (3, 0.)],
                limits(),
                electric_controls(),
                &mut || false,
            )
            .unwrap();
        near(continuing.field().voltage_v()[1], 3. * a * a);
        near(
            continuing.field().stored_energy_j() + continuing.step().heat_j + next.step().heat_j,
            9.,
        );
        near(old.voltage_v()[1], 3.);
        assert_eq!(old.time_s(), 0.);
    }
}

#[test]
fn native_interface_charge_and_electrode_work_follow_heterogeneous_materials() {
    let left = material(0., 1., false);
    let right = material(4., 1., false);
    let context = point(310., 0.5);
    let sources = [
        conductor(&left, &context, false),
        conductor(&right, &context, false),
    ];
    let mesh = TriComplex2::from_indexed_triangles(
        "two leaky dielectric layers",
        vec![[0., 0.], [1., 0.], [2., 0.], [0., 1.], [1., 1.], [2., 1.]],
        vec![[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4]],
        Metric2::planar(1.).unwrap(),
    )
    .unwrap();
    let conduction = PlanarConduction::new(
        mesh,
        vec![
            sources[0].clone(),
            sources[0].clone(),
            sources[1].clone(),
            sources[1].clone(),
        ],
        limits(),
        &mut || false,
    )
    .unwrap();
    let old = PlanarLeakyDielectric::new(
        conduction,
        vec![0., 2., 4., 0., 2., 4.],
        &[0, 3, 2, 5],
        0.,
        limits(),
        electric_controls(),
        &mut || false,
    )
    .unwrap();
    let next = old
        .step_to_candidate(
            1.,
            old.conduction().conductors().to_vec(),
            &[(0, 0.), (3, 0.), (2, 4.), (5, 4.)],
            limits(),
            electric_controls(),
            &mut || false,
        )
        .unwrap();
    // Each slab has epsilon=2 F/m, but sigma=[1,1/3] S/m. The interface
    // voltage obeys 4*u' + (4/3)*u = 4/3. Starting u=2, midpoint gives u1=12/7.
    // Interface charge q=4*u-8 therefore changes from zero to -8/7 C.
    near(next.field().voltage_v()[1], 12. / 7.);
    near(next.field().voltage_v()[4], 12. / 7.);
    near(old.nodal_charge_c()[1] + old.nodal_charge_c()[4], 0.);
    near(
        next.field().nodal_charge_c()[1] + next.field().nodal_charge_c()[4],
        -8. / 7.,
    );
    near(next.step().stored_before_j, 8.);
    near(next.step().stored_after_j, 400. / 49.);
    near(next.step().heat_j, 244. / 49.);
    near(next.step().supplied_work_j, 36. / 7.);
    let terminal = next.step().nodal_boundary_current_a[2].unwrap()
        + next.step().nodal_boundary_current_a[5].unwrap();
    near(terminal, 9. / 7.);
    for i in 0..4 {
        near(
            next.step().cell_power_w[i],
            if i < 2 { 169. / 98. } else { 75. / 98. },
        );
    }
    // Nonzero interface storage and unequal currents are actual spatial
    // redistribution; replacing the whole body by a rescaled R cannot supply it.
    near(next.step().conduction_current_density_a_m2[0][0], -13. / 7.);
    near(next.step().conduction_current_density_a_m2[2][0], -5. / 7.);
    assert!(next.step().energy_residual_j.abs() < electric_controls().energy_tolerance_j);
    let ramped = next
        .field()
        .step_to_candidate(
            2.,
            next.field().conduction().conductors().to_vec(),
            &[(0, 0.), (3, 0.), (2, 5.), (5, 5.)],
            limits(),
            electric_controls(),
            &mut || false,
        )
        .unwrap();
    // The same interface equation now includes epsilon_right*dV/dt. An
    // endpoint-only resistive source omits this 2 A displacement contribution.
    let u0 = 12. / 7.;
    let u1 = 387. / 196.;
    let middle = 0.5 * (u0 + u1);
    let current = 2. * (1. - (u1 - u0)) + (4.5 - middle) / 3.;
    near(ramped.field().voltage_v()[1], u1);
    near(
        ramped.step().heat_j,
        middle * middle + (4.5 - middle).powi(2) / 3.,
    );
    near(ramped.step().supplied_work_j, 4.5 * current);
}

fn thermal_reference(v0: f64, h0: f64, slope: f64) -> (f64, f64) {
    // Unit volume, epsilon=2, dt=1. Exact storage loss h1-h0=v0²-v1²
    // supplies latent enthalpy. Midpoint relaxation uses R(f1)=1+slope*f1.
    let (mut low, mut high) = (0., v0);
    for _ in 0..90 {
        let v = 0.5 * (low + high);
        let h = h0 + v0 * v0 - v * v;
        let resistance = 1. + slope * (h - 20.) / 100.;
        let rhs = v0 * (4. * resistance - 1.) / (4. * resistance + 1.);
        if v < rhs {
            low = v;
        } else {
            high = v;
        }
    }
    let v = 0.5 * (low + high);
    (v, h0 + v0 * v0 - v * v)
}

#[test]
fn stored_spatial_charge_melts_material_and_phase_data_changes_its_relaxation() {
    let mut remaining = Vec::new();
    for tensor in [false, true] {
        for slope in [0., 3.] {
            let material = material(slope, 1., false);
            let caloric = curve(&material, 1.);
            let thermal = initial(&caloric, 20.);
            let context = point(310., 0.);
            let old = field(&material, &context, tensor, 8.);
            let cells = [binding(&material, &caloric, &context)];
            let layout = network(&cells);
            let first = layout
                .step_dielectric_candidate(
                    &thermal,
                    1.,
                    &old,
                    &[(0, 0.), (3, 0.)],
                    limits(),
                    electric_controls(),
                    config(),
                    &[25.],
                    || false,
                )
                .unwrap();
            let (v1, h1) = thermal_reference(8., 20., slope);
            near(first.field().voltage_v()[1], v1);
            near(first.thermal().energy_j(), h1);
            near(first.thermal().temperature_k(CellId(7)).unwrap(), 310.);
            near(
                first.report().source.cells[0].phase.liquid_mass_fraction(),
                (h1 - 20.) / 100.,
            );
            near(first.report().deposited_joule_j, 64. - v1 * v1);
            near(first.report().electrical.supplied_work_j, 0.);
            near(
                first.thermal().energy_j() + first.field().stored_energy_j(),
                84.,
            );
            check_sources(
                old.conduction(),
                first.field().conduction(),
                310.,
                (h1 - 20.) / 100.,
                slope,
            );
            let second = layout
                .step_dielectric_candidate(
                    first.thermal(),
                    2.,
                    first.field(),
                    &[(0, 0.), (3, 0.)],
                    limits(),
                    electric_controls(),
                    config(),
                    &[10.],
                    || false,
                )
                .unwrap();
            let (v2, h2) = thermal_reference(v1, h1, slope);
            near(second.field().voltage_v()[1], v2);
            near(second.thermal().energy_j(), h2);
            near(
                second.thermal().energy_j() + second.field().stored_energy_j(),
                84.,
            );
            assert_eq!(second.field().time_s(), second.thermal().time_s());
            for source in second.field().dielectric().materials() {
                near(source.query().axes()[FRACTION], (h2 - 20.) / 100.);
                let epsilon = source.property(ELECTRIC_PERMITTIVITY_PROPERTY).unwrap();
                near(epsilon.value_si(), 2.);
                let original = old.dielectric().materials()[0]
                    .property(ELECTRIC_PERMITTIVITY_PROPERTY)
                    .unwrap();
                assert_eq!(
                    epsilon.answer().receipt.selected,
                    original.answer().receipt.selected
                );
            }
            near(thermal.energy_j(), 20.);
            near(old.stored_energy_j(), 64.);
            remaining.push(v2);
        }
    }
    assert!(remaining[1] > remaining[0]);
    near(remaining[0], remaining[2]);
    near(remaining[1], remaining[3]);
}

#[test]
fn changing_dielectric_storage_missing_source_and_late_refusal_publish_no_state() {
    let context = point(310., 0.);
    let changing = material(3., 1., true);
    let conduction = PlanarConduction::new(
        mesh(),
        vec![conductor(&changing, &context, false); 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(
        PlanarLeakyDielectric::new(
            conduction,
            vec![0., 8., 8., 0., 4.],
            &[0, 3],
            0.,
            limits(),
            electric_controls(),
            &mut || false
        )
        .is_err()
    );
    let missing = card(3., 1.);
    assert!(
        PlanarLeakyDielectric::new(
            problem(&missing, &context, false),
            vec![0., 8., 8., 0., 4.],
            &[0, 3],
            0.,
            limits(),
            electric_controls(),
            &mut || false
        )
        .is_err()
    );
    let material = material(3., 1., false);
    let caloric = curve(&material, 1.);
    let thermal = initial(&caloric, 20.);
    let old = field(&material, &context, true, 8.);
    let original_charge = old.nodal_charge_c().to_vec();
    let cells = [binding(&material, &caloric, &context)];
    let layout = network(&cells);
    let mut polls = 0;
    let accepted = layout
        .step_dielectric_candidate(
            &thermal,
            1.,
            &old,
            &[(0, 0.), (3, 0.)],
            limits(),
            electric_controls(),
            config(),
            &[25.],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(accepted.thermal().energy_j() > thermal.energy_j());
    let mut calls = 0;
    assert!(matches!(
        layout.step_dielectric_candidate(
            &thermal,
            1.,
            &old,
            &[(0, 0.), (3, 0.)],
            limits(),
            electric_controls(),
            config(),
            &[25.],
            || {
                calls += 1;
                calls >= polls
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert!(
        old.step_to_candidate(
            1.,
            old.conduction().conductors().to_vec(),
            &[(0, 0.), (3, 0.), (1, 8.), (2, 8.)],
            limits(),
            electric_controls(),
            &mut || false
        )
        .is_err()
    );
    let restricted = self::material(3., 0.1, false);
    let restricted_curve = curve(&restricted, 1.);
    let restricted_thermal = initial(&restricted_curve, 20.);
    let restricted_field = field(&restricted, &context, true, 8.);
    let restricted_cell = binding(&restricted, &restricted_curve, &context);
    assert!(
        network(&[restricted_cell])
            .step_dielectric_candidate(
                &restricted_thermal,
                1.,
                &restricted_field,
                &[(0, 0.), (3, 0.)],
                limits(),
                electric_controls(),
                config(),
                &[5.],
                || false
            )
            .is_err()
    );
    assert_eq!(old.nodal_charge_c(), original_charge.as_slice());
    assert_eq!(old.time_s(), 0.);
    assert_eq!(thermal.time_s(), 0.);
    near(old.stored_energy_j(), 64.);
    near(thermal.energy_j(), 20.);
}
