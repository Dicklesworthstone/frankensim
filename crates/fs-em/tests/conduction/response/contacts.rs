//! G1/G3: nonlinear bulk laws and native contact jumps share one spatial field.
use super::*;
use fs_em::conduction::ConductionContact;
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
    resolve_electrical_interface_state,
};

fn contact(
    a: &ScalarConductorResponse,
    b: &ScalarConductorResponse,
    resistance: f64,
) -> ConductionContact {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
                ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Scalar {
                value: resistance,
                dims: ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured constant contact".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: a.original().material().material().clone(),
            texture_frame: "left".into(),
        },
        SurfaceSpec {
            material: b.original().material().material().clone(),
            texture_frame: "right".into(),
        },
        SystemContext {
            medium: "manufactured".into(),
            third_body: None,
            environment: "fixed".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    ConductionContact {
        first_vertices: [1, 2],
        second_vertices: [4, 7],
        source: resolve_electrical_interface_state(
            &card,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    }
}

#[test]
fn nonlinear_series_field_retains_surface_current_heat_and_source_identity() {
    let a = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        1.0,
        0.2,
        10.0,
    );
    let b = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        3.0,
        0.4,
        10.0,
    );
    let mut currents = Vec::new();
    for rc in [0.2, 1.3] {
        let mesh = TriComplex2::from_indexed_triangles(
            "nonlinear contact specimen",
            vec![
                [0.0, 0.0],
                [1.0, 0.0],
                [1.0, 1.0],
                [0.0, 1.0],
                [1.0, 0.0],
                [2.0, 0.0],
                [2.0, 1.0],
                [1.0, 1.0],
            ],
            vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
            Metric2::planar(1.0).unwrap(),
        )
        .unwrap();
        let model = PlanarConductionResponse::new_with_phase_contacts(
            mesh,
            vec![a.clone(), a.clone(), b.clone(), b.clone()],
            vec![contact(&a, &b, rc)],
            vec![[None, None]],
            limits(),
            &mut || false,
        )
        .unwrap();
        let fields = |j: f64| {
            [
                2.0 * j / (1.0 + (1.0 + 0.8 * j).sqrt()),
                2.0 * j / (3.0 + (9.0 + 1.6 * j).sqrt()),
            ]
        };
        let (mut lo, mut hi) = (0.0, 10.0);
        for _ in 0..100 {
            let j = 0.5 * (lo + hi);
            if fields(j).iter().sum::<f64>() + rc * j < 2.0 {
                lo = j;
            } else {
                hi = j;
            }
        }
        let j = 0.5 * (lo + hi);
        let e = fields(j);
        currents.push(j);
        for sign in [1.0, -1.0] {
            let voltage = [(0, 0.0), (3, 0.0), (5, sign * 2.0), (6, sign * 2.0)];
            let solved = model
                .solve_candidate(&voltage, &[], limits(), controls(), policy(), &mut || false)
                .unwrap();
            let report = solved.report();
            assert!(report.field_solves > 1);
            for cell in 0..4 {
                near(report.current_density_a_m2[cell][0], -sign * j);
                near(report.electric_field_v_m[cell][0], -sign * e[cell / 2]);
                near(report.cell_power_w[cell], 0.5 * j * e[cell / 2]);
            }
            near(report.contact_response[0].current_a, -sign * j);
            near(report.contact_response[0].voltage_jump_v[0], -sign * j * rc);
            near(report.contact_response[0].power_w, j * j * rc);
            near(report.power_w, 2.0 * j);
            near(
                report.cell_power_w.iter().sum::<f64>() + report.contact_response[0].power_w,
                report.power_w,
            );
            let current = model
                .solve_current_candidate(
                    &[5, 6],
                    &[0, 3],
                    sign * j,
                    limits(),
                    controls(),
                    policy(),
                    &mut || false,
                )
                .unwrap();
            near(current.report().voltage_v[5], sign * 2.0);
            let (native, sources) = solved.into_parts().0.into_frozen_parts();
            let continued =
                PlanarConductionResponse::from_frozen_parts(native, sources, &mut || false)
                    .unwrap();
            assert_eq!(continued.contacts(), model.contacts());
            assert!(
                continued
                    .solve_with_materials_and_contacts_candidate(
                        continued.conductors().to_vec(),
                        vec![contact(&a, &b, rc + 1.0)],
                        &voltage,
                        &[],
                        limits(),
                        controls(),
                        policy(),
                        &mut || false,
                    )
                    .is_err()
            );
        }
    }
    assert!(currents[0] > currents[1]);
}
