//! General layered response preserves the original constitutive and spatial owners.
use super::*;
use fs_material::visco_tensor::kelvin_control::TensorKelvinStressControlBudget;

fn viscous_elastic() -> TensorProny {
    let v = [0.2, -0.1, 0.3, 0.15, -0.05, 0.4];
    let d = core::array::from_fn(|i| core::array::from_fn(|j| {
        v[i] * v[j] + if i == j { 0.05 } else { 0.0 }
    }));
    TensorProny::new_mandel(elastic(0.03).stiffness(), vec![], 0.03)
        .unwrap().with_kelvin_viscosity_mandel(d).unwrap()
        .with_band((0.0, 20.0)).unwrap()
}

fn mixed_section(kelvin: TensorProny) -> LaminateSection {
    LaminateSection::try_new(vec![
        layer(LayerLaw::Prony(kelvin), 0.08, 4.0, 0.31),
        layer(LayerLaw::Prony(prony()), 0.02, 9.0, -0.47),
    ], ReferenceSurface::Offset(0.03)).unwrap()
}

#[test]
fn general_engine_preserves_maxwell_only_physical_results_and_channels() {
    with_cx(|cx| {
        let shell = source(laminate(), CstDktMass::Consistent, cx).unwrap();
        let legacy = shell.prepare_material_causal_observer(causal_budget(), cx).unwrap();
        let general = shell.prepare_viscoelastic_material_causal_observer(
            causal_budget(), TensorKelvinStressControlBudget::default(), cx,
        ).unwrap();
        let old = legacy.start_relaxed(Dd::ZERO, cx).unwrap();
        let new = general.start_relaxed(Dd::ZERO, cx).unwrap();
        let q = motion();
        let banks = q.each_ref().map(|q| q.as_slice());
        let a = legacy.advance_window(&old, &FREQUENCIES, &banks, &[0.0; 6], 0.2, cx).unwrap();
        let b = general.advance_window(&new, &FREQUENCIES, &banks, &[0.0; 6], 0.2, cx).unwrap();
        assert_eq!(a.history.full_displacement(), b.history.full_displacement());
        assert_eq!(a.history.time_s_dd(), b.history.time_s_dd());
        assert_eq!(a.field.material_heat_j, b.field.material_heat_j);
        assert_eq!(a.field.stress_work_j, b.field.stress_work_j);
        assert_eq!(a.field.stored_energy_change_j, b.field.stored_energy_change_j);
        assert_eq!(a.field.kelvin_heat_j, b.field.kelvin_heat_j.unwrap());
        assert_eq!(a.field.maxwell_heat_j, b.field.maxwell_heat_j.unwrap());
        for (a, b) in a.field.points.iter().zip(&b.field.points) {
            assert_eq!(a.end_strain, b.end_strain);
            assert_eq!(a.end_stress, b.end_stress);
            assert_eq!(a.strain_bounds, b.strain_bounds);
            assert_eq!(a.material.material_heat_j_m3, b.material.material_heat_j_m3);
            assert_eq!(b.end_strain_rate, None);
        }
        assert!(general.advance_window(&old, &FREQUENCIES, &banks, &[0.0; 6], 0.2, cx).is_err());
    });
}

#[test]
fn kelvin_and_maxwell_layers_reconstruct_physical_stress_and_retain_heat_on_hold() {
    with_cx(|cx| {
        let kelvin = viscous_elastic();
        let section = mixed_section(kelvin.clone());
        let shell = source(section.clone(), CstDktMass::Consistent, cx).unwrap();
        assert!(shell.prepare_material_causal_observer(causal_budget(), cx).is_err());
        let observer = shell.prepare_viscoelastic_material_causal_observer(
            causal_budget(), TensorKelvinStressControlBudget::default(), cx,
        ).unwrap();
        let initial = observer.start_relaxed(Dd::ZERO, cx).unwrap();
        let q = motion();
        let banks = q.each_ref().map(|q| q.as_slice());
        let step = observer.advance_window(&initial, &FREQUENCIES, &banks, &[0.0; 6], 0.2, cx).unwrap();
        assert_eq!(step.field.kelvin_heat_j, None);
        assert_eq!(step.field.maxwell_heat_j, None);
        assert!(step.field.material_heat_j > 0.0);
        let mut deposited = 0.0;
        let mut layer_heat = [0.0; 2];
        let weights = [1.0, 1.0, 1.0, std::f64::consts::SQRT_2,
            std::f64::consts::SQRT_2, std::f64::consts::SQRT_2];
        for point in &step.field.points {
            let layer = &section.layers()[point.layer];
            for k in 0..2 {
                let strain = principal(&q[k], point.triangle_barycentric,
                    point.thickness_position_m, layer.orientation_rad);
                for j in 0..3 { near(point.in_plane_strain[k][j], strain[j], 3e-12); }
            }
            deposited += point.volume_m3 * point.material.material_heat_j_m3;
            layer_heat[point.layer] += point.volume_m3 * point.material.material_heat_j_m3;
            if point.layer == 0 {
                let rate = point.end_strain_rate.unwrap();
                // No Maxwell branch in this layer: reconstruct every physical
                // stress directly from the ORIGINAL full C and D tensors.
                for i in 0..6 {
                    let expected: f64 = (0..6).map(|j| {
                        (kelvin.c_eq()[i][j] * point.end_strain[j]
                            + kelvin.kv_pa_s().unwrap()[i][j] * rate[j]) * weights[j] / weights[i]
                    }).sum();
                    assert!((point.end_stress[i] - expected).abs() < 2e-12);
                }
                assert!(point.material.kelvin_heat_j_m3.is_none());
            } else {
                assert_eq!(point.material.kelvin_heat_j_m3, Some(0.0));
                assert!(point.material.maxwell_heat_j_m3.unwrap() > 0.0);
            }
            for i in [2, 4, 5] { assert_eq!(point.end_stress[i], 0.0); }
        }
        assert!(layer_heat.iter().all(|q| *q > 0.0));
        real_near(step.field.material_heat_j, deposited, 3e-13);
        let balance = step.field.stress_work_j - step.field.stored_energy_change_j
            - step.field.material_heat_j;
        assert!(balance.abs() < 3e-10 * step.field.stress_work_j.abs());
        let held = observer.advance_window(&step.history, &[], &[],
            &step.history.full_displacement()[12..], 0.125, cx).unwrap();
        assert_eq!(held.field.stress_work_j, 0.0);
        assert!(held.field.material_heat_j > 0.0);
        real_near(held.field.stored_energy_change_j, -held.field.material_heat_j, 3e-10);
        assert_eq!(held.history.full_displacement(), step.history.full_displacement());
        assert_eq!(initial.stored_energy_j(), 0.0);
        assert_eq!(initial.full_displacement(), &[0.0; 18]);
    });
}

#[test]
fn kelvin_controller_capacity_refuses_but_original_directional_viscosity_evolves() {
    with_cx(|cx| {
        let shell = source(mixed_section(viscous_elastic()), CstDktMass::Consistent, cx).unwrap();
        let too_small = TensorKelvinStressControlBudget {
            max_state_coordinates: 2,
            ..TensorKelvinStressControlBudget::default()
        };
        assert!(shell.prepare_viscoelastic_material_causal_observer(causal_budget(), too_small, cx).is_err());
        let mut d = [[0.0; 6]; 6];
        d[0][0] = 0.1;
        let singular = prony().with_kelvin_viscosity_mandel(d).unwrap();
        let shell = source(mixed_section(singular), CstDktMass::Consistent, cx).unwrap();
        // The original in-plane dashpot has exactly zero transverse viscosity.
        // Its elastic traction closure and every Maxwell memory remain physical.
        let observer = shell.prepare_viscoelastic_material_causal_observer(
            causal_budget(), TensorKelvinStressControlBudget::default(), cx,
        ).unwrap();
        let initial = observer.start_relaxed(Dd::ZERO, cx).unwrap();
        let q = motion();
        let banks = q.each_ref().map(|q| q.as_slice());
        let step = observer.advance_window(
            &initial, &FREQUENCIES, &banks, &[0.0; 6], 0.2, cx,
        ).unwrap();
        assert!(step.field.material_heat_j > 0.0);
        assert!(step.field.points.iter().any(|point| {
            point.layer == 0 && point.material.material_heat_j_m3 > 0.0
        }));
        for point in &step.field.points {
            for i in [2, 4, 5] { assert_eq!(point.end_stress[i], 0.0); }
        }
        let scale = step.field.stress_work_j.abs()
            + step.field.stored_energy_change_j.abs() + step.field.material_heat_j;
        assert!(step.field.energy_balance_residual_j.abs() <= 3e-10 * scale);
        assert_eq!(initial.time_s_dd(), Dd::ZERO);
        assert_eq!(initial.dissipated_energy_j(), 0.0);
    });
}
