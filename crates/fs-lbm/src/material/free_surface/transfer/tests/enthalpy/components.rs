//! Detached liquid pieces retain their own motion and heat through native ingress.
use super::*;

const CENTRES: [[f64; 3]; 2] = [[0.004; 3], [0.012, 0.004, 0.004]];
const TRANSLATIONS: [[f64; 3]; 2] = [[0.01, -0.005, 0.002], [-0.007, 0.004, -0.003]];

fn pieces() -> (Vec<[f64; 3]>, Vec<[u32; 4]>, Vec<[f64; 3]>, Vec<f64>) {
    let (mut p, mut t, _) = cube([0.0; 3]);
    let (second, topology, _) = cube([0.008, 0.0, 0.0]);
    p.extend(second);
    t.extend(topology.into_iter().map(|tet| tet.map(|node| node + 8)));
    let v = p
        .iter()
        .enumerate()
        .map(|(index, &point)| {
            let piece = index / 8;
            add(
                TRANSLATIONS[piece],
                cross(
                    [0.0, 0.0, if piece == 0 { 3.0 } else { -3.0 }],
                    sub(point, CENTRES[piece]),
                ),
            )
        })
        .collect();
    let h = (0..16)
        .map(|index| if index < 8 { 150000.0 } else { 300000.0 })
        .collect();
    (p, t, v, h)
}

fn bounds() -> TetEnthalpyTransferLimits {
    let mut result = thermal_limits();
    result.mechanical.max_tetrahedra = 16;
    result
}

#[test]
fn g1_counter_rotating_heated_fragments_preserve_independent_momenta_and_continue_native_flow() {
    let card = thermal_card(1000.0, 2e-5);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point(),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let surface = resolve_liquid_surface_state(
        &card,
        &point(),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let caloric = chart(&card, 1000.0);
    let (p, t, v, h) = pieces();
    let unchanged = (p.clone(), t.clone(), v.clone(), h.clone());
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, 0.001, 0.001, 0.1).unwrap();
        let (mut thermal, report) = scaling
            .free_surface3_from_tetrahedra_with_enthalpy(
                domain(),
                [0.0; 3],
                source(&p, &t, &v, &h),
                &law,
                &caloric,
                bounds(),
                FreeSurfaceEnthalpyControls::default(),
                cx,
                &gate,
            )
            .unwrap();
        assert_eq!(report.mechanical.components.len(), 2);
        assert_eq!(report.components.len(), 2);
        let mass = 1000.0 * 0.004_f64.powi(3);
        let inertia = mass * 0.004_f64.powi(2) / 6.0;
        let discrete_inertia = inertia * 15.0 / 16.0;
        let expected_change = 0.5 * 9.0 * inertia * (inertia / discrete_inertia - 1.0);
        for piece in 0..2 {
            let component = &report.mechanical.components[piece];
            let p = scale(TRANSLATIONS[piece], mass);
            let omega = if piece == 0 { 3.0 } else { -3.0 };
            let l = add(cross(CENTRES[piece], p), [0.0, 0.0, inertia * omega]);
            assert_eq!(
                component.source_tetrahedra,
                (piece * 6..piece * 6 + 6).collect::<Vec<_>>()
            );
            assert_eq!(component.target_cells.len(), 64);
            close(component.source.mass_kg, mass, 2e-18);
            close(component.target.mass_kg, mass, 2e-18);
            close(component.kinetic_energy_change_j, expected_change, 4e-22);
            for axis in 0..3 {
                close(component.source.momentum_kg_m_s[axis], p[axis], 3e-20);
                close(component.target.momentum_kg_m_s[axis], p[axis], 3e-20);
                close(
                    component.source.angular_momentum_kg_m2_s[axis],
                    l[axis],
                    5e-22,
                );
                close(
                    component.target.angular_momentum_kg_m2_s[axis],
                    l[axis],
                    5e-22,
                );
                close(component.centroid_shift_m[axis], 0.0, 2e-16);
            }
            let expected_h = mass * if piece == 0 { 150000.0 } else { 300000.0 };
            close(
                report.components[piece].source_enthalpy_j,
                expected_h,
                2e-13,
            );
            close(
                report.components[piece].target_enthalpy_j,
                expected_h,
                2e-12,
            );
            close(report.components[piece].enthalpy_defect_j, 0.0, 2e-12);
            for &at in &component.target_cells {
                assert_eq!(
                    thermal.material_at(at).unwrap().material().card_identity(),
                    card.content_hash()
                );
            }
        }
        // Each component agrees with its independent one-body transfer. The
        // remote counter-rotating piece cannot supply or absorb its torque.
        let (single, _) = scaling
            .free_surface3_from_tetrahedra_with_enthalpy(
                domain(),
                [0.0; 3],
                source(&p[..8], &t[..6], &v[..8], &h[..8]),
                &law,
                &caloric,
                bounds(),
                FreeSurfaceEnthalpyControls::default(),
                cx,
                &gate,
            )
            .unwrap();
        for &[x, y, z] in &report.mechanical.components[0].target_cells {
            let alone = single.flow().sample_at_si(x, y, z).unwrap().unwrap();
            let together = thermal.flow().sample_at_si(x, y, z).unwrap().unwrap();
            for axis in 0..3 {
                close(together.velocity_m_s[axis], alone.velocity_m_s[axis], 2e-15);
            }
        }
        let step = thermal.step(&gate).unwrap();
        close(step.mass_after_kg, 2.0 * mass, 1e-14);
        close(step.enthalpy_after_j, mass * 450000.0, 1e-9);
        close(step.external_heat_j, 0.0, 0.0);
        close(step.viscous_heat_j, 0.0, 0.0);
        assert_eq!(thermal.accepted_steps(), 1);
        close(thermal.physical_time_s(), 0.001, 1e-16);
    });
    assert_eq!((p, t, v, h), unchanged);
}

#[test]
fn g3_opposite_component_projection_energy_errors_cannot_cancel_the_budget() {
    let material = liquid(1000.0);
    let (p, t, mut v, _) = pieces();
    let strain_rate = 3.0 * (16.0_f64 / 15.0).sqrt();
    for (index, velocity) in v.iter_mut().enumerate() {
        let piece = index / 8;
        let r = sub(p[index], CENTRES[piece]);
        *velocity = if piece == 0 {
            cross([0.0, 0.0, 3.0], r)
        } else {
            [strain_rate * r[0], -strain_rate * r[1], 0.0]
        };
    }
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&material, 0.001, 0.001, 0.1).unwrap();
        let input = TetrahedralLiquidState {
            positions_m: &p,
            tetrahedra: &t,
            velocities_m_s: &v,
        };
        let (_, report) = scaling
            .free_surface3_from_tetrahedra(domain(), [0.0; 3], input, bounds().mechanical, cx)
            .unwrap();
        let mass = 1000.0 * 0.004_f64.powi(3);
        let expected = mass * strain_rate.powi(2) * 0.001_f64.powi(2) / 12.0;
        close(
            report.components[0].kinetic_energy_change_j,
            expected,
            2e-22,
        );
        close(
            report.components[1].kinetic_energy_change_j,
            -expected,
            2e-22,
        );
        close(report.kinetic_energy_change_j, 0.0, 4e-22);
        let mut tight = bounds().mechanical;
        tight.maximum_projection_energy_change_j = 1.5 * expected;
        assert!(matches!(
            scaling.free_surface3_from_tetrahedra(domain(), [0.0; 3], input, tight, cx,),
            Err(TetLiquidTransferError::Invalid(
                "component transfer exceeds the shared absolute conservation/error budget"
            ))
        ));
    });
}

#[test]
fn g4_unresolved_disconnected_pieces_refuse_without_merging_or_publishing_heat() {
    let card = thermal_card(1000.0, 2e-5);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point(),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let surface = resolve_liquid_surface_state(
        &card,
        &point(),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let caloric = chart(&card, 1000.0);
    let mut p = Vec::new();
    for shift in [0.0021, 0.0026] {
        p.extend([
            [shift, 0.0021, 0.0021],
            [shift + 0.0001, 0.0021, 0.0021],
            [shift, 0.0022, 0.0021],
            [shift, 0.0021, 0.0022],
        ]);
    }
    let t = [[0, 1, 2, 3], [4, 5, 6, 7]];
    let v = [[0.0; 3]; 8];
    let h = [
        100000.0, 100000.0, 100000.0, 100000.0, 300000.0, 300000.0, 300000.0, 300000.0,
    ];
    let before = (p.clone(), t, v, h);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, 0.001, 0.001, 0.1).unwrap();
        assert!(matches!(
            scaling.free_surface3_from_tetrahedra_with_enthalpy(
                domain(),
                [0.0; 3],
                source(&p, &t, &v, &h),
                &law,
                &caloric,
                bounds(),
                FreeSurfaceEnthalpyControls::default(),
                cx,
                &gate,
            ),
            Err(TetEnthalpyTransferError::Transfer(
                TetLiquidTransferError::Invalid(
                    "disconnected source pieces share a receiving cell; refine before transfer"
                )
            ))
        ));
    });
    assert_eq!((p, t, v, h), before);
    let (p, t, v, h) = pieces();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    with_cx(&cancelled, |cx| {
        let scaling = plan_free_surface_from_material(&surface, 0.001, 0.001, 0.1).unwrap();
        assert!(
            scaling
                .free_surface3_from_tetrahedra_with_enthalpy(
                    domain(),
                    [0.0; 3],
                    source(&p, &t, &v, &h),
                    &law,
                    &caloric,
                    bounds(),
                    FreeSurfaceEnthalpyControls::default(),
                    cx,
                    &cancelled,
                )
                .is_err()
        );
    });
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, 0.001, 0.001, 0.1).unwrap();
        assert!(
            scaling
                .free_surface3_from_tetrahedra_with_enthalpy(
                    domain(),
                    [0.0; 3],
                    source(&p, &t, &v, &h),
                    &law,
                    &caloric,
                    bounds(),
                    FreeSurfaceEnthalpyControls::default(),
                    cx,
                    &gate,
                )
                .is_ok()
        );
    });
}
