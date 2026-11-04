use super::*;
use fs_material::thermovisco::{Limits, ThermalKnot};

fn law(branches: Vec<(f64, f64)>) -> ThermoMaxwell {
    ThermoMaxwell::new(1e8, branches, [300.0, 500.0].into_iter().map(|t| ThermalKnot {
        temperature_k: t, modulus_scale: 1.0, relaxation_scale: 1.0, free_strain: (t-300.0)*1e-5,
    }).collect(), Limits { max_strain: 0.03, max_strain_rate_s: 10.0 }).unwrap()
}
fn strip<'a>() -> LayeredStrip<'a> {
    LayeredStrip::new(0.01, Budget { max_layers: 32, max_branch_states: 1024 }).unwrap()
}
fn geometry(z: f64, thickness: f64) -> LayerGeometry {
    LayerGeometry { area_m2: 0.002*thickness, centroid_z_m: z, thickness_m: thickness }
}
fn config() -> Config {
    Config { force_tolerance_n: 1e-8, moment_tolerance_nm: 1e-10, energy_tolerance_j: 1e-9,
        max_strain: 0.03, max_rotation_rad: 0.1 }
}
fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() < tol, "{a:e} != {b:e}"); }

#[test]
fn uniform_free_contraction_and_fixed_reaction_match_analytic_bar() {
    let law = law(vec![]);
    let strip = strip().deposit(1, geometry(0.0, 0.001), &law, 500.0).unwrap();
    let t = BTreeMap::from([(1, 300.0)]);
    let (free, r) = strip.advance_to(1.0, &t, Control::RELEASED, config(), || Ok(())).unwrap();
    close(free.axial_strain(), -0.002, 1e-14); close(free.curvature_per_m(), 0.0, 1e-12);
    close(free.deformed_position_m(0.01, 0.0).unwrap()[0], 0.00998, 1e-14);
    close(r.axial_force_n, 0.0, 1e-8);
    let (_, held) = strip.advance_to(1.0, &t, Control::HELD, config(), || Ok(())).unwrap();
    close(held.axial_reaction_n, 0.4, 1e-12); close(held.bending_reaction_nm, 0.0, 1e-12);
    close(held.energy_j.mechanical_work, 0.0, 1e-12);
    for report in [r, held] {
        assert_eq!(
            report.axial_force_n.to_bits(),
            report.work_conjugate_axial_force_n.to_bits()
        );
        assert_eq!(
            report.bending_moment_nm.to_bits(),
            report.work_conjugate_bending_moment_nm.to_bits()
        );
    }
}

#[test]
fn differential_layer_cooling_bends_after_prescribed_release() {
    let law = law(vec![]);
    let initial = strip().deposit(1, geometry(-0.0005, 0.001), &law, 500.0).unwrap()
        .deposit(2, geometry(0.0005, 0.001), &law, 500.0).unwrap();
    let t = BTreeMap::from([(1, 300.0), (2, 400.0)]);
    let (held, reactions) = initial.advance_to(1.0, &t, Control::HELD, config(), || Ok(())).unwrap();
    assert!(reactions.axial_reaction_n > 0.0 && reactions.bending_reaction_nm > 0.0);
    assert_eq!(held.deformed_position_m(0.01, 0.0).unwrap()[1], 0.0);
    let (released, result) = held.advance_to(1.1, &t, Control::RELEASED, config(), || Ok(())).unwrap();
    close(released.axial_strain(), -0.0015, 1e-14);
    close(released.curvature_per_m(), -0.75, 1e-12);
    close(released.deformed_position_m(0.01, 0.0).unwrap()[1], -3.75e-5, 1e-14);
    assert!(result.energy_j.stored_change < 0.0);
    assert!(result.force_residual_n.abs() < 1e-8 && result.moment_residual_nm.abs() < 1e-10);
    assert!(result.energy_j.residual.abs() < 1e-9);
}

#[test]
fn new_layer_uses_current_deformed_reference_without_resetting_old_stress() {
    let law = law(vec![(3e8, 2.0)]);
    let strip = strip().deposit(1, geometry(0.0, 0.001), &law, 400.0).unwrap();
    let control = Control { axial: AxialControl::Strain(0.001), bending: BendingControl::CurvaturePerM(0.2) };
    let (loaded, _) = strip.advance_to(1.0, &BTreeMap::from([(1, 350.0)]), control, config(), || Ok(())).unwrap();
    let born = loaded.deposit(2, geometry(0.001, 0.001), &law, 450.0).unwrap();
    assert_eq!(loaded.layers()[&1], born.layers()[&1]);
    for (point, z) in born.layers()[&2].material_points().iter().zip(born.layers()[&2].geometry.points()) {
        assert_eq!(point.stress_pa(), 0.0);
        close(point.birth_strain(), loaded.axial_strain()-z*loaded.curvature_per_m(), 1e-15);
    }
    assert_eq!(born.layers()[&2].birth_time_s, 1.0);
}

#[test]
fn rectangle_refinement_preserves_axial_bending_and_history_energy() {
    let law = law(vec![(3e8, 2.0)]);
    let a = strip().deposit(1, geometry(0.0, 0.002), &law, 400.0).unwrap();
    let b = strip().deposit(1, geometry(-0.0005, 0.001), &law, 400.0).unwrap()
        .deposit(2, geometry(0.0005, 0.001), &law, 400.0).unwrap();
    let load = Control { axial: AxialControl::Strain(0.001), bending: BendingControl::CurvaturePerM(0.2) };
    let (a, ar) = a.advance_to(1.0, &BTreeMap::from([(1, 350.0)]), load, config(), || Ok(())).unwrap();
    let (b, br) = b.advance_to(1.0, &BTreeMap::from([(1, 350.0), (2, 350.0)]), load, config(), || Ok(())).unwrap();
    close(ar.axial_force_n, br.axial_force_n, 1e-12); close(ar.bending_moment_nm, br.bending_moment_nm, 1e-12);
    close(a.stored_energy_j(), b.stored_energy_j(), 1e-12);
    close(ar.energy_j.dissipated, br.energy_j.dissipated, 1e-12);
}

#[test]
fn cold_relaxation_keeps_force_balance_and_irreversible_history_after_release() {
    let law = law(vec![(3e8, 2.0)]);
    let initial = strip().deposit(1, geometry(0.0, 0.001), &law, 500.0).unwrap();
    let t = BTreeMap::from([(1, 300.0)]);
    let (held, a) = initial.advance_to(1.0, &t, Control::HELD, config(), || Ok(())).unwrap();
    let (relaxed, b) = held.advance_to(3.0, &t, Control::HELD, config(), || Ok(())).unwrap();
    assert!(b.axial_reaction_n < a.axial_reaction_n && b.energy_j.dissipated > 0.0);
    let (free, r) = relaxed.advance_to(3.1, &t, Control::RELEASED, config(), || Ok(())).unwrap();
    assert!(free.axial_strain() < 0.0 && free.axial_strain() > -0.002);
    assert!(free.layers()[&1].material_points()[0].branch_strains()[0].abs() > 1e-6);
    close(r.axial_force_n, 0.0, 1e-8);
}

#[test]
fn cancellation_and_large_rotation_leave_original_checkpoint_unchanged() {
    let law = law(vec![(3e8, 2.0)]);
    let initial = strip().deposit(1, geometry(0.0, 0.001), &law, 400.0).unwrap();
    let old = initial.clone(); let t = BTreeMap::from([(1, 350.0)]);
    let mut checkpoints = 0;
    assert!(initial.advance_to(1.0, &t, Control::HELD, config(), || {
        checkpoints += 1; if checkpoints >= 3 { Err(Error::Cancelled) } else { Ok(()) }
    }).is_err());
    assert_eq!(initial, old);
    let excessive = Control { axial: AxialControl::Strain(0.0), bending: BendingControl::CurvaturePerM(20.0) };
    assert!(matches!(initial.advance_to(1.0, &t, excessive, config(), || Ok(())), Err(Error::Regime(_))));
    assert_eq!(initial, old);
}

#[test]
fn exact_elastic_strain_ramp_uses_interval_work_without_numerical_heat() {
    let law = law(vec![])
        .with_history_integration(fs_material::thermovisco::HistoryIntegration::ExactStrainRamp);
    let geometry = geometry(0.0, 0.002);
    let original = strip().deposit(1, geometry, &law, 300.0).unwrap();
    let control = Control {
        axial: AxialControl::Strain(0.002),
        bending: BendingControl::CurvaturePerM(0.4),
    };
    let (next, report) = original
        .advance_to(
            1.0,
            &BTreeMap::from([(1, 300.0)]),
            control,
            config(),
            || Ok(()),
        )
        .unwrap();
    let ea = 1e8 * geometry.area_m2;
    let ei = ea * geometry.thickness_m.powi(2) / 12.0;
    let expected = 0.5 * original.length_m() * (ea * 0.002_f64.powi(2) + ei * 0.4_f64.powi(2));
    close(report.axial_force_n, ea * 0.002, 1e-12);
    close(report.bending_moment_nm, ei * 0.4, 1e-14);
    close(report.work_conjugate_axial_force_n, 0.5 * ea * 0.002, 1e-12);
    close(
        report.work_conjugate_bending_moment_nm,
        0.5 * ei * 0.4,
        1e-14,
    );
    close(next.stored_energy_j(), expected, 1e-14);
    close(report.energy_j.mechanical_work, expected, 1e-14);
    assert_eq!(report.energy_j.dissipated, 0.0);
    assert_eq!(report.energy_j.numerical_dissipation, 0.0);
    assert_eq!(original.axial_strain(), 0.0);
}

#[test]
fn exact_maxwell_work_and_birth_preserve_the_original_relaxation_history() {
    let law = law(vec![(3e8, 2.0)]).with_history_integration(HistoryIntegration::ExactStrainRamp);
    let shape = geometry(0.0, 0.002);
    let original = strip().deposit(1, shape, &law, 300.0).unwrap();
    let control = Control {
        axial: AxialControl::Strain(0.002),
        bending: BendingControl::CurvaturePerM(0.4),
    };
    let (loaded, report, layers) = original
        .advance_with_layer_energy(
            0.4,
            &BTreeMap::from([(1, 300.0)]),
            control,
            config(),
            || Ok(()),
        )
        .unwrap();
    // Independent Maxwell solution q(t)=strain_rate*tau*(1-exp(-t/tau)).
    let ratio: f64 = 0.4 / 2.0;
    let endpoint_factor = -(-ratio).exp_m1() / ratio;
    let mean_factor = (1.0 - endpoint_factor) / ratio;
    let strain_square_volume = original.length_m()
        * shape.area_m2
        * (0.002_f64.powi(2) + shape.thickness_m.powi(2) * 0.4_f64.powi(2) / 12.0);
    let expected_work = (0.5e8 + 3e8 * mean_factor) * strain_square_volume;
    let expected_storage = 0.5 * (1e8 + 3e8 * endpoint_factor.powi(2)) * strain_square_volume;
    close(report.energy_j.mechanical_work, expected_work, 1e-13);
    close(loaded.stored_energy_j(), expected_storage, 1e-13);
    close(
        layers[&1].dissipated,
        expected_work - expected_storage,
        1e-13,
    );
    assert_eq!(report.energy_j.numerical_dissipation, 0.0);
    let newborn = loaded
        .deposit(2, geometry(0.0015, 0.001), &law, 300.0)
        .unwrap();
    assert_eq!(newborn.layers()[&1], loaded.layers()[&1]);
    assert_eq!(newborn.stored_energy_j(), loaded.stored_energy_j());
    assert!(
        newborn.layers()[&2]
            .material_points()
            .iter()
            .all(|p| p.stress_pa() == 0.0)
    );
    let (relaxed, hold) = newborn
        .advance_to(
            2.4,
            &BTreeMap::from([(1, 300.0), (2, 300.0)]),
            control,
            config(),
            || Ok(()),
        )
        .unwrap();
    let expected_heat =
        0.5 * 3e8 * endpoint_factor.powi(2) * strain_square_volume * (1.0 - (-2.0_f64).exp());
    close(hold.energy_j.dissipated, expected_heat, 1e-13);
    close(
        newborn.stored_energy_j() - relaxed.stored_energy_j(),
        expected_heat,
        1e-13,
    );
    assert_eq!(hold.energy_j.mechanical_work, 0.0);
    assert_eq!(hold.energy_j.numerical_dissipation, 0.0);
}

#[test]
fn explicit_mean_force_and_moment_use_the_original_exact_ramp_tangent() {
    let law = law(vec![(3e8, 2.0)]).with_history_integration(HistoryIntegration::ExactStrainRamp);
    let geometry = geometry(0.0, 0.002);
    let original = strip().deposit(1, geometry, &law, 300.0).unwrap();
    let dt = 0.4;
    let ratio: f64 = dt / 2.0;
    let mean_modulus = 0.5e8 + 3e8 * (1.0 + (-ratio).exp_m1() / ratio) / ratio;
    let force = 0.1;
    let moment = 1e-5;
    let expected_strain = force / (mean_modulus * geometry.area_m2);
    let expected_curvature =
        moment / (mean_modulus * geometry.area_m2 * geometry.thickness_m.powi(2) / 12.0);
    for axial in [
        AxialControl::MeanForceN(force),
        AxialControl::Strain(expected_strain),
    ] {
        let (next, report) = original
            .advance_to(
                dt,
                &BTreeMap::from([(1, 300.0)]),
                Control {
                    axial,
                    bending: BendingControl::MeanMomentNm(moment),
                },
                config(),
                || Ok(()),
            )
            .unwrap();
        close(next.axial_strain(), expected_strain, 1e-14);
        close(next.curvature_per_m(), expected_curvature, 1e-12);
        close(report.work_conjugate_axial_force_n, force, 1e-12);
        close(report.work_conjugate_bending_moment_nm, moment, 1e-14);
        assert!(report.axial_force_n > force && report.bending_moment_nm > moment);
        close(
            report.energy_j.mechanical_work,
            original.length_m() * (force * expected_strain + moment * expected_curvature),
            1e-13,
        );
        assert_eq!(report.energy_j.numerical_dissipation, 0.0);
        assert!(report.energy_j.residual.abs() < config().energy_tolerance_j);
    }
}

#[test]
fn incompatible_load_histories_refuse_without_changing_deposited_material() {
    let exact = law(vec![(3e8, 2.0)]).with_history_integration(HistoryIntegration::ExactStrainRamp);
    let backward = law(vec![(3e8, 2.0)]);
    let original = strip()
        .deposit(1, geometry(0.0, 0.001), &exact, 300.0)
        .unwrap();
    let saved = original.clone();
    let temperatures = BTreeMap::from([(1, 300.0)]);
    for control in [
        Control::RELEASED,
        Control {
            axial: AxialControl::Strain(0.0),
            bending: BendingControl::MomentNm(0.0),
        },
        Control {
            axial: AxialControl::ForceN(0.0),
            bending: BendingControl::MeanMomentNm(0.0),
        },
    ] {
        assert!(matches!(
            original.advance_to(1.0, &temperatures, control, config(), || Ok(())),
            Err(Error::Invalid(_))
        ));
    }
    assert_eq!(original, saved);
    let mixed = original
        .deposit(2, geometry(0.001, 0.001), &backward, 300.0)
        .unwrap();
    assert!(matches!(
        mixed.advance_to(
            1.0,
            &BTreeMap::from([(1, 300.0), (2, 300.0)]),
            Control {
                axial: AxialControl::MeanForceN(0.1),
                bending: BendingControl::CurvaturePerM(0.0)
            },
            config(),
            || Ok(())
        ),
        Err(Error::Invalid(_))
    ));
}
