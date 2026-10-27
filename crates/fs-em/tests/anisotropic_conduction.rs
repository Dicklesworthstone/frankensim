//! Independent field/current/power checks for source-oriented electrical data.
mod common;
use common::{config, ends, layered, limits, near};
use fs_em::conduction::{ConductionControls, ConductionCurrentTrace, PlanarConduction};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::conductor::planar::{
    PRINCIPAL_CONDUCTIVITY_PROPERTIES, PRINCIPAL_RESISTIVITY_PROPERTIES, PlanarElectricalFrame,
    ResolvedOrthotropicConductor2, ResolvedPlanarConductor, resolve_orthotropic_conductor2,
};
use fs_material::conductor::{ElectricalPropertyBasis, resolve_conductor_state_with_basis};
use fs_material::state_point::MaterialPropertySelection;
use fs_rep_mesh::{Metric2, TriComplex2};

fn card(names: &[&str], values: &[f64], basis: ElectricalPropertyBasis) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (&name, &value) in names.iter().zip(values) {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, basis.dimensions()),
                value: PropertyValue::Scalar {
                    value,
                    dims: basis.dimensions(),
                },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic directional conductor".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic electrical tensor".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn frame(angle: f64) -> PlanarElectricalFrame {
    // Distinct declared fixture frame identities, independent from specimen data.
    let a = card(
        &PRINCIPAL_CONDUCTIVITY_PROPERTIES,
        &[17., 19.],
        ElectricalPropertyBasis::Conductivity,
    );
    let b = card(
        &PRINCIPAL_CONDUCTIVITY_PROPERTIES,
        &[23., 29.],
        ElectricalPropertyBasis::Conductivity,
    );
    PlanarElectricalFrame::new(a.content_hash(), b.content_hash(), angle).unwrap()
}
fn conductor(
    sigma: [f64; 2],
    basis: ElectricalPropertyBasis,
    angle: f64,
) -> ResolvedOrthotropicConductor2 {
    let names = match basis {
        ElectricalPropertyBasis::Conductivity => PRINCIPAL_CONDUCTIVITY_PROPERTIES,
        ElectricalPropertyBasis::Resistivity => PRINCIPAL_RESISTIVITY_PROPERTIES,
    };
    let values = sigma.map(|v| match basis {
        ElectricalPropertyBasis::Conductivity => v,
        ElectricalPropertyBasis::Resistivity => v.recip(),
    });
    resolve_orthotropic_conductor2(
        &card(&names, &values, basis),
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        basis,
        frame(angle),
    )
    .unwrap()
}
fn square(depth: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "tensor-conductor-square",
        vec![[0., 0.], [1., 0.], [1., 1.], [0., 1.], [0.5, 0.5]],
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap()
}
fn controls() -> ConductionControls {
    ConductionControls {
        linear: config(),
        current_tolerance_a: 1e-10,
        power_tolerance_w: 1e-10,
    }
}

#[test]
fn oriented_sources_produce_transverse_current_and_independent_joule_work() {
    for basis in [
        ElectricalPropertyBasis::Conductivity,
        ElectricalPropertyBasis::Resistivity,
    ] {
        let source = conductor([2., 8.], basis, core::f64::consts::FRAC_PI_4);
        let sigma = source.conductivity_s_m();
        near(sigma[0][0], 5.);
        near(sigma[1][1], 5.);
        near(sigma[0][1], -3.);
        let problem =
            PlanarConduction::new(square(0.3), vec![source; 4], limits(), &mut || false).unwrap();
        // V=x+2y has zero divergence of constant tensor current, including at
        // the unprescribed interior vertex. E=(-1,-2), J=(1,-7), E.J=13.
        for (sign, offset) in [(1., 0.), (-1., 7.)] {
            let field = problem
                .solve(
                    &[
                        (0, offset),
                        (1, offset + sign),
                        (2, offset + 3. * sign),
                        (3, offset + 2. * sign),
                    ],
                    controls(),
                    &mut || false,
                )
                .unwrap();
            near(field.voltage_v()[4], offset + 1.5 * sign);
            for i in 0..4 {
                near(field.electric_field_v_m()[i][0], -sign);
                near(field.electric_field_v_m()[i][1], -2. * sign);
                near(field.current_density_a_m2()[i][0], sign);
                near(field.current_density_a_m2()[i][1], -7. * sign);
                near(field.cell_power_w()[i], 13. * 0.3 / 4.);
            }
            near(field.power_w(), 13. * 0.3);
            near(field.boundary_current_a(), 0.);
            near(field.power_balance_residual_w(), 0.);
        }
    }
}

#[test]
fn tensor_current_traces_solve_the_actual_free_voltages() {
    let source = conductor(
        [2., 8.],
        ElectricalPropertyBasis::Conductivity,
        core::f64::consts::FRAC_PI_4,
    );
    let problem =
        PlanarConduction::new(square(0.4), vec![source; 4], limits(), &mut || false).unwrap();
    // V=x+y, J=(-2,-2). Inward loads are +2 on right/top and -2
    // on bottom. Only the full left trace has prescribed potential.
    let currents = [
        ConductionCurrentTrace {
            vertices: [1, 2],
            inward_current_density_a_m2: 2.,
        },
        ConductionCurrentTrace {
            vertices: [2, 3],
            inward_current_density_a_m2: 2.,
        },
        ConductionCurrentTrace {
            vertices: [0, 1],
            inward_current_density_a_m2: -2.,
        },
    ];
    let field = problem
        .solve_with_current_traces(&[(0, 0.), (3, 1.)], &currents, controls(), &mut || false)
        .unwrap();
    for (&actual, expected) in field.voltage_v().iter().zip([0., 1., 2., 1., 1.]) {
        near(actual, expected);
    }
    for j in field.current_density_a_m2() {
        near(j[0], -2.);
        near(j[1], -2.);
    }
    near(field.power_w(), 1.6);
    near(field.power_balance_residual_w(), 0.);
}

#[test]
fn mixed_scalar_and_oriented_cells_share_one_resistance_and_heat_operator() {
    let scalar_card = card(
        &["electrical_conductivity"],
        &[4.],
        ElectricalPropertyBasis::Conductivity,
    );
    let scalar = resolve_conductor_state_with_basis(
        &scalar_card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        ElectricalPropertyBasis::Conductivity,
    )
    .unwrap();
    for (angle, sigma_x) in [(0., 2.), (core::f64::consts::FRAC_PI_2, 8.)] {
        let anisotropic = conductor([2., 8.], ElectricalPropertyBasis::Resistivity, angle);
        let materials: Vec<ResolvedPlanarConductor> = vec![
            anisotropic.clone().into(),
            anisotropic.into(),
            scalar.clone().into(),
            scalar.clone().into(),
        ];
        let problem =
            PlanarConduction::new(layered(0.5, 1., 0.2), materials, limits(), &mut || false)
                .unwrap();
        let expected_r = (0.5 / sigma_x + 0.5 / 4.) / 0.2;
        let resistance = problem
            .resistance(&[2, 5], &[0, 3], controls(), &mut || false)
            .unwrap();
        near(resistance.resistance_ohm(), expected_r);
        let field = resistance
            .field_at_voltage(3., 1e-10, 1e-10, &mut || false)
            .unwrap();
        near(field.power_w(), 9. / expected_r);
        for (i, j) in field.current_density_a_m2().iter().enumerate() {
            near(j[0], -3. / expected_r / 0.2);
            near(j[1], 0.);
            let sigma = if i < 2 { sigma_x } else { 4. };
            near(field.cell_power_w()[i], j[0] * j[0] / sigma * 0.05);
        }
    }
    // Isotropic degeneration still enters the same solver with tensor data.
    let isotropic = conductor([4., 4.], ElectricalPropertyBasis::Conductivity, 0.71);
    let tensor = PlanarConduction::new(
        layered(0.5, 1., 0.2),
        vec![isotropic; 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let native = PlanarConduction::new(
        layered(0.5, 1., 0.2),
        vec![scalar; 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let a = tensor.solve(&ends(2.), controls(), &mut || false).unwrap();
    let b = native.solve(&ends(2.), controls(), &mut || false).unwrap();
    for (a, b) in a
        .current_density_a_m2()
        .iter()
        .flatten()
        .zip(b.current_density_a_m2().iter().flatten())
    {
        near(*a, *b);
    }
    near(a.power_w(), b.power_w());
}

#[test]
fn incomplete_axes_and_inconsistent_geometry_frames_refuse() {
    let missing = card(
        &PRINCIPAL_CONDUCTIVITY_PROPERTIES[..1],
        &[2.],
        ElectricalPropertyBasis::Conductivity,
    );
    assert!(
        resolve_orthotropic_conductor2(
            &missing,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            ElectricalPropertyBasis::Conductivity,
            frame(0.)
        )
        .is_err()
    );
    assert!(
        PlanarElectricalFrame::new(missing.content_hash(), missing.content_hash(), 0.2).is_err()
    );
    let a = conductor([2., 8.], ElectricalPropertyBasis::Conductivity, 0.);
    let bad_frame =
        PlanarElectricalFrame::new(frame(0.).material_frame(), missing.content_hash(), 0.).unwrap();
    let b = ResolvedOrthotropicConductor2::try_from_material(
        a.material().clone(),
        a.source_basis(),
        bad_frame,
    )
    .unwrap();
    assert!(
        PlanarConduction::new(
            square(1.),
            vec![a.clone(), a.clone(), b, a],
            limits(),
            &mut || false
        )
        .is_err()
    );
}
