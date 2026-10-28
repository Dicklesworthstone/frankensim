//! Independent motional-EMF, anisotropic current, heat and work oracles.
mod common;
use common::{config, ends, layered, limits, near};
use fs_em::conduction::motion::ConductionMotion;
use fs_em::conduction::{
    ConductionControls, ConductionCurrentTrace, ConductionError, PlanarConduction,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::conductor::planar::{
    PRINCIPAL_CONDUCTIVITY_PROPERTIES, PRINCIPAL_RESISTIVITY_PROPERTIES, PlanarElectricalFrame,
    ResolvedOrthotropicConductor2, resolve_orthotropic_conductor2,
};
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
fn open_translating_anisotropic_bar_screens_motional_emf_without_braking() {
    for basis in [
        ElectricalPropertyBasis::Conductivity,
        ElectricalPropertyBasis::Resistivity,
    ] {
        for (angle, b, offset) in [(0.0, 0.5, 0.0), (0.7, -0.5, 7.0)] {
            let source = conductor([2.0, 8.0], basis, angle);
            let problem =
                PlanarConduction::new(square(0.3), vec![source; 4], limits(), &mut || false)
                    .unwrap();
            let velocity = [[2.0, 0.0]; 4];
            let field = problem
                .solve_with_motion(
                    &[(0, offset), (1, offset)],
                    &[],
                    ConductionMotion {
                        velocity_m_s: &velocity,
                        normal_flux_density_t: &[b; 4],
                    },
                    controls(),
                    &mut || false,
                )
                .unwrap();
            // Only the bottom reference trace is held. All other exterior
            // edges insulate; the solver must recover V=offset-u B y.
            for (&voltage, position) in field.voltage_v().iter().zip(problem.mesh().vertices()) {
                near(voltage, offset - 2.0 * b * position[1]);
            }
            for cell in 0..4 {
                near(field.electric_field_v_m()[cell][0], 0.0);
                near(field.electric_field_v_m()[cell][1], 2.0 * b);
                for axis in 0..2 {
                    near(field.conductor_field_v_m()[cell][axis], 0.0);
                    near(field.current_density_a_m2()[cell][axis], 0.0);
                    near(field.lorentz_force_density_n_m3()[cell][axis], 0.0);
                }
            }
            assert_eq!(field.velocity_m_s(), velocity);
            assert_eq!(field.normal_flux_density_t(), &[b; 4]);
            near(field.joule_power_w(), 0.0);
            near(field.mechanical_power_w(), 0.0);
            near(field.terminal_power_w(), 0.0);
            near(field.boundary_current_a(), 0.0);
        }
    }
}

#[test]
fn shorted_loaded_generator_and_motor_obey_directional_conductivity_and_work() {
    // Rotating the same material by 90 degrees exchanges the axial values.
    for (angle, sigma_y) in [(0.0, 8.0), (core::f64::consts::FRAC_PI_2, 2.0)] {
        let source = conductor([2.0, 8.0], ElectricalPropertyBasis::Conductivity, angle);
        let problem =
            PlanarConduction::new(square(0.25), vec![source; 4], limits(), &mut || false).unwrap();
        for e_y in [0.0, 0.4, 1.8] {
            let field = problem
                .solve_with_motion(
                    &[(0, 0.0), (1, 0.0), (2, -e_y), (3, -e_y)],
                    &[],
                    ConductionMotion {
                        velocity_m_s: &[[3.0, 0.0]; 4],
                        normal_flux_density_t: &[0.4; 4],
                    },
                    controls(),
                    &mut || false,
                )
                .unwrap();
            // Independent uniform-bar formula. E_y=0 shorts the electrodes;
            // 0<E_y<uB generates into a load; E_y>uB drives the moving bar.
            let effective_e_y = e_y - 1.2;
            let j_y = sigma_y * effective_e_y;
            let heat = 0.25 * sigma_y * effective_e_y * effective_e_y;
            let mechanical = 0.25 * 1.2 * j_y;
            let electrical = 0.25 * e_y * j_y;
            near(field.voltage_v()[4], -0.5 * e_y);
            for cell in 0..4 {
                near(field.current_density_a_m2()[cell][0], 0.0);
                near(field.current_density_a_m2()[cell][1], j_y);
                near(field.lorentz_force_density_n_m3()[cell][0], 0.4 * j_y);
                near(field.lorentz_force_density_n_m3()[cell][1], 0.0);
                near(field.cell_joule_power_w()[cell], heat / 4.0);
                near(field.cell_mechanical_power_w()[cell], mechanical / 4.0);
            }
            near(field.joule_power_w(), heat);
            near(field.mechanical_power_w(), mechanical);
            near(field.electrical_power_w(), electrical);
            near(field.terminal_power_w(), electrical);
            near(
                field.joule_power_w() + field.mechanical_power_w(),
                field.terminal_power_w(),
            );
            near(field.maximum_free_current_residual_a(), 0.0);
            near(field.power_balance_residual_w(), 0.0);
            near(field.work_pairing_residual_w(), 0.0);
            if e_y == 0.0 {
                assert!(field.mechanical_power_w() < 0.0);
                near(field.joule_power_w(), -field.mechanical_power_w());
            } else if e_y < 1.2 {
                assert!(field.terminal_power_w() < 0.0);
            } else {
                assert!(field.mechanical_power_w() > 0.0);
            }
        }
    }
}

#[test]
fn rotated_tensor_pairs_transverse_force_with_conductor_frame_heat() {
    let source = conductor(
        [2.0, 8.0],
        ElectricalPropertyBasis::Resistivity,
        core::f64::consts::FRAC_PI_4,
    );
    let problem =
        PlanarConduction::new(square(0.3), vec![source; 4], limits(), &mut || false).unwrap();
    let field = problem
        .solve_with_motion(
            &[(0, 0.0), (1, 1.0), (2, 3.0), (3, 2.0)],
            &[],
            ConductionMotion {
                velocity_m_s: &[[0.5, 1.0]; 4],
                normal_flux_density_t: &[2.0; 4],
            },
            controls(),
            &mut || false,
        )
        .unwrap();
    // Sigma=[[5,-3],[-3,5]], E=(-1,-2), u x B=(2,-1).
    // Therefore E'=(1,-3), J=(14,-18), f=(-36,-28),
    // q=68, mechanical work density=-46, electrical work density=22.
    near(field.voltage_v()[4], 1.5);
    for cell in 0..4 {
        for (axis, expected) in [14.0, -18.0].into_iter().enumerate() {
            near(field.current_density_a_m2()[cell][axis], expected);
        }
        near(field.conductor_field_v_m()[cell][0], 1.0);
        near(field.conductor_field_v_m()[cell][1], -3.0);
        near(field.lorentz_force_density_n_m3()[cell][0], -36.0);
        near(field.lorentz_force_density_n_m3()[cell][1], -28.0);
        near(field.cell_joule_power_w()[cell], 68.0 * 0.3 / 4.0);
        near(field.cell_mechanical_power_w()[cell], -46.0 * 0.3 / 4.0);
    }
    near(field.joule_power_w(), 68.0 * 0.3);
    near(field.mechanical_power_w(), -46.0 * 0.3);
    near(field.terminal_power_w(), 22.0 * 0.3);
    near(field.electrical_power_w(), 22.0 * 0.3);
    near(field.boundary_current_a(), 0.0);
    near(field.power_balance_residual_w(), 0.0);
}

#[test]
fn imposed_current_retains_heat_and_shifts_terminal_voltage_by_motional_emf() {
    let source = conductor([2.0, 8.0], ElectricalPropertyBasis::Conductivity, 0.0);
    let problem =
        PlanarConduction::new(square(0.4), vec![source; 4], limits(), &mut || false).unwrap();
    let currents = [ConductionCurrentTrace {
        vertices: [2, 3],
        inward_current_density_a_m2: 2.0,
    }];
    let boundary = [(0, 5.0), (1, 5.0)];
    let stationary = problem
        .solve_with_current_traces(&boundary, &currents, controls(), &mut || false)
        .unwrap();
    let field = problem
        .solve_with_motion(
            &boundary,
            &currents,
            ConductionMotion {
                velocity_m_s: &[[2.0, 0.0]; 4],
                normal_flux_density_t: &[0.5; 4],
            },
            controls(),
            &mut || false,
        )
        .unwrap();
    // The current source fixes J_y=-2 A/m2. The motional emf changes the
    // terminal voltage by -u B L=-1 V, but leaves J.Sigma^-1.J unchanged.
    near(field.voltage_v()[2], 4.25);
    near(field.voltage_v()[3], 4.25);
    near(field.voltage_v()[2] - stationary.voltage_v()[2], -1.0);
    for j in field.current_density_a_m2() {
        near(j[0], 0.0);
        near(j[1], -2.0);
    }
    near(field.joule_power_w(), 0.2);
    near(field.joule_power_w(), stationary.power_w());
    near(field.mechanical_power_w(), -0.8);
    near(field.terminal_power_w(), -0.6);
    near(field.boundary_current_a(), 0.0);
    near(field.maximum_free_current_residual_a(), 0.0);
    assert_eq!(
        field.nodal_applied_current_a(),
        stationary.nodal_applied_current_a()
    );
}

#[test]
fn heterogeneous_countermoving_layers_share_one_continuous_current() {
    let left = conductor([2.0, 7.0], ElectricalPropertyBasis::Conductivity, 0.0);
    let right = conductor(
        [3.0, 8.0],
        ElectricalPropertyBasis::Resistivity,
        core::f64::consts::FRAC_PI_2,
    );
    let problem = PlanarConduction::new(
        layered(0.5, 0.6, 0.4),
        vec![left.clone(), left, right.clone(), right],
        limits(),
        &mut || false,
    )
    .unwrap();
    let field = problem
        .solve_with_motion(
            &ends(0.0),
            &[],
            ConductionMotion {
                velocity_m_s: &[[0.0, 1.0], [0.0, 1.0], [0.0, -2.0], [0.0, -2.0]],
                normal_flux_density_t: &[0.5; 4],
            },
            controls(),
            &mut || false,
        )
        .unwrap();
    // Shorted series layers: J=(L1*m1+L2*m2)/(L1/sigma1+L2/sigma2)
    // =(.5*.5-.5*1)/(.5/2+.5/8)=-.8. Enforcing current locally as
    // Sigma*(u x B) would violate continuity and give the wrong force/heat.
    near(field.voltage_v()[1], 0.45);
    near(field.voltage_v()[4], 0.45);
    for cell in 0..4 {
        near(field.current_density_a_m2()[cell][0], -0.8);
        near(field.current_density_a_m2()[cell][1], 0.0);
        near(
            field.electric_field_v_m()[cell][0],
            if cell < 2 { -0.9 } else { 0.9 },
        );
        near(field.lorentz_force_density_n_m3()[cell][1], 0.4);
        near(
            field.cell_joule_power_w()[cell],
            if cell < 2 { 0.0192 } else { 0.0048 },
        );
        near(
            field.cell_mechanical_power_w()[cell],
            if cell < 2 { 0.024 } else { -0.048 },
        );
    }
    near(field.joule_power_w(), 0.048);
    near(field.mechanical_power_w(), -0.048);
    near(field.terminal_power_w(), 0.0);
    near(field.work_pairing_residual_w(), 0.0);
    near(field.power_balance_residual_w(), 0.0);
}

#[test]
fn zero_motional_emf_preserves_stationary_solve_and_stationary_lorentz_force() {
    let source = conductor(
        [2.0, 8.0],
        ElectricalPropertyBasis::Conductivity,
        core::f64::consts::FRAC_PI_4,
    );
    let problem =
        PlanarConduction::new(square(0.3), vec![source; 4], limits(), &mut || false).unwrap();
    let boundary = [(0, 0.0), (1, 1.0), (2, 3.0), (3, 2.0)];
    let stationary = problem.solve(&boundary, controls(), &mut || false).unwrap();
    for (velocity, b) in [([0.5, 1.0], 0.0), ([0.0, 0.0], 0.4)] {
        let field = problem
            .solve_with_motion(
                &boundary,
                &[],
                ConductionMotion {
                    velocity_m_s: &[velocity; 4],
                    normal_flux_density_t: &[b; 4],
                },
                controls(),
                &mut || false,
            )
            .unwrap();
        assert_eq!(field.voltage_v(), stationary.voltage_v());
        assert_eq!(field.electric_field_v_m(), stationary.electric_field_v_m());
        assert_eq!(
            field.current_density_a_m2(),
            stationary.current_density_a_m2()
        );
        assert_eq!(field.cell_joule_power_w(), stationary.cell_power_w());
        assert_eq!(
            field.nodal_boundary_current_a(),
            stationary.nodal_boundary_current_a()
        );
        near(field.joule_power_w(), stationary.power_w());
        near(field.mechanical_power_w(), 0.0);
        for force in field.lorentz_force_density_n_m3() {
            near(force[0], -7.0 * b);
            near(force[1], -b);
        }
    }
    for motion in [
        ConductionMotion {
            velocity_m_s: &[[0.0, 0.0]; 3],
            normal_flux_density_t: &[0.0; 4],
        },
        ConductionMotion {
            velocity_m_s: &[[f64::NAN, 0.0]; 4],
            normal_flux_density_t: &[0.0; 4],
        },
        ConductionMotion {
            velocity_m_s: &[[0.0, 0.0]; 4],
            normal_flux_density_t: &[f64::INFINITY; 4],
        },
    ] {
        assert!(matches!(
            problem.solve_with_motion(&boundary, &[], motion, controls(), &mut || false),
            Err(ConductionError::Invalid(_))
        ));
    }
    assert!(matches!(
        problem.solve_with_motion(
            &boundary,
            &[],
            ConductionMotion {
                velocity_m_s: &[[0.0, 0.0]; 4],
                normal_flux_density_t: &[0.4; 4]
            },
            controls(),
            &mut || true,
        ),
        Err(ConductionError::Cancelled)
    ));
}
