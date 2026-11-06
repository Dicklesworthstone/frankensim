//! Native implicit SGS convection keeps one physical heat and momentum clock.
use super::*;
use fs_lbm::conjugate::{
    EddyWallTreatment, FvSgsConfig, TransientConfig, march_conjugate_with_sgs_implicit_buoyancy,
    march_energy,
};
use fs_lbm::sgs3::{SgsModel, SgsUse};

fn sgs(model: SgsModel) -> FvSgsConfig {
    FvSgsConfig {
        model,
        filter_width_factor: 1.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

#[test]
fn g1_each_sgs_model_closes_buoyancy_and_uses_its_projected_heat_field_once() {
    let (domain, mut flow, thermal, _, initial) = cavity();
    flow.turbulent_prandtl = 0.6;
    let solids = [SolidMaterial::new("nonlinear-wall", 1.0)
        .with_heat_capacity(2.0)
        .with_conductivity_table(&[(300.0, 1.0), (330.0, 2.5)])];
    let b = buoyancy(0.001, [0.0, -9.81, 0.0]);
    let clock = clock(0.1, 1);
    let fluid = fluid();
    for model in [
        SgsModel::Smagorinsky { cs: 0.5 },
        SgsModel::Wale { cw: 0.5 },
        SgsModel::Vreman { cv: 0.5 },
    ] {
        let result = march_conjugate_with_sgs_implicit_buoyancy(
            &domain,
            &fluid,
            &solids,
            &thermal,
            &flow,
            &clock,
            sgs(model),
            &b,
            controls(),
            &initial,
            |_| 1.0,
            |_| 1.0,
            &energy(),
            None,
            &CancelGate::new(),
        )
        .unwrap();
        let march = &result.solution;
        let projected = &march.flow.flow;
        assert!(projected.eddy_viscosity_m2_s.iter().any(|&nu| nu > 0.0));
        assert!(
            projected
                .heat_eddy_viscosity_m2_s
                .iter()
                .any(|&nu| nu > 0.0)
        );
        assert!(
            projected
                .velocity_m_s
                .iter()
                .flatten()
                .any(|u| u.abs() > 1e-8)
        );
        assert_eq!(projected.turbulent_prandtl, 0.6);
        assert!(result.records[0].passes > 1);
        assert!(result.records[0].temperature_defect_k <= controls().temperature_tolerance_k);
        assert!(result.records[0].force_defect_n_m3 <= controls().force_tolerance_n_m3);
        assert_eq!(march.flow.records.len(), 1);
        assert_eq!(march.energy_records.len(), 1);
        assert_eq!(march.flow.averaged_steps, 1);

        // Recover the same single physical heat interval independently from
        // its accepted projected flow. Form k_t explicitly so a hidden Pr_t
        // fallback or accumulation over nonlinear passes cannot satisfy this.
        let mut explicit = thermal.clone();
        explicit.eddy_wall_treatment = EddyWallTreatment::Molecular;
        explicit.eddy_conductivity_w_m_k = projected
            .heat_eddy_viscosity_m2_s
            .iter()
            .map(|nu_t| fluid.density_kg_m3 * fluid.specific_heat_j_kg_k * nu_t / 0.6)
            .collect();
        for (cell, &conductivity) in explicit.eddy_conductivity_w_m_k.iter().enumerate() {
            if !domain.is_fluid(cell) {
                assert_eq!(conductivity, 0.0);
            }
        }
        let heat_clock = TransientConfig {
            time_step_s: clock.time_step_s,
            steps: 1,
            energy: energy(),
        };
        let frozen = march_energy(
            &domain,
            &fluid,
            &solids,
            &projected.field,
            &explicit,
            &initial,
            |_| 1.0,
            &heat_clock,
            None,
            &CancelGate::new(),
        )
        .unwrap();
        let molecular = march_energy(
            &domain,
            &fluid,
            &solids,
            &projected.field,
            &thermal,
            &initial,
            |_| 1.0,
            &heat_clock,
            None,
            &CancelGate::new(),
        )
        .unwrap();
        let difference = |other: &[f64]| {
            march
                .temperature
                .iter()
                .zip(other)
                .fold(0.0_f64, |largest, (a, b)| largest.max((a - b).abs()))
        };
        assert!(
            difference(&frozen.temperature) < 1e-8,
            "SGS heat must use its current projected field, including native solid k(T)"
        );
        assert!(
            difference(&molecular.temperature) > 1e-10,
            "the selected SGS model must affect convection-diffusion"
        );
        let stored: f64 = initial
            .iter()
            .zip(&march.temperature)
            .enumerate()
            .map(|(cell, (&old, &new))| 2.0 * domain.volume(cell) * (new - old))
            .sum();
        assert!(stored.abs() < 1e-8);
        assert!(march.energy_records[0].closure_j.abs() < 1e-8);
        assert_eq!(march.energy_records[0].source_j, 0.0);
        assert_eq!(march.energy_records[0].boundary_outflow_j, 0.0);
    }
}

#[test]
fn g3_zero_sgs_is_the_exact_native_implicit_march_with_the_same_history() {
    let (domain, flow, thermal, solids, initial) = cavity();
    let b = buoyancy(0.001, [0.0, -9.81, 0.0]);
    let clock = clock(0.1, 3);
    let laminar = march_conjugate_implicit_buoyancy(
        &domain,
        &fluid(),
        &solids,
        &thermal,
        &flow,
        &clock,
        &b,
        controls(),
        &initial,
        |_| 1.0,
        |_| 1.0,
        &energy(),
        None,
        &CancelGate::new(),
    )
    .unwrap();
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let disabled = march_conjugate_with_sgs_implicit_buoyancy(
            &domain,
            &fluid(),
            &solids,
            &thermal,
            &flow,
            &clock,
            sgs(model),
            &b,
            controls(),
            &initial,
            |_| 1.0,
            |_| 1.0,
            &energy(),
            None,
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(disabled.solution, laminar.solution);
        assert_eq!(disabled.records, laminar.records);
        assert_eq!(disabled.momentum_solves, laminar.momentum_solves);
        assert_eq!(disabled.energy_solves, laminar.energy_solves);
        assert_eq!(disabled.source_refreshes, laminar.source_refreshes);
    }
}

#[test]
fn g4_native_implicit_sgs_retains_model_cap_and_competing_closure_refusals() {
    let (domain, flow, thermal, solids, initial) = cavity();
    let b = buoyancy(0.001, [0.0, -9.81, 0.0]);
    let clock = clock(0.1, 1);
    let model = sgs(SgsModel::Smagorinsky { cs: 0.5 });
    for rejected in [
        FvSgsConfig {
            usage: SgsUse::AcceptanceRender,
            ..model
        },
        FvSgsConfig {
            maximum_eddy_to_molecular: 0.0,
            ..model
        },
    ] {
        assert!(matches!(
            march_conjugate_with_sgs_implicit_buoyancy(
                &domain,
                &fluid(),
                &solids,
                &thermal,
                &flow,
                &clock,
                rejected,
                &b,
                controls(),
                &initial,
                |_| 1.0,
                |_| 1.0,
                &energy(),
                None,
                &CancelGate::new(),
            ),
            Err(ChtError::InvalidInput {
                field: "simple.sgs",
                ..
            })
        ));
    }
    let mut thermal = thermal.clone();
    thermal.eddy_conductivity_w_m_k = vec![0.1; domain.cell_count()];
    assert!(matches!(
        march_conjugate_with_sgs_implicit_buoyancy(
            &domain,
            &fluid(),
            &solids,
            &thermal,
            &flow,
            &clock,
            model,
            &b,
            controls(),
            &initial,
            |_| 1.0,
            |_| 1.0,
            &energy(),
            None,
            &CancelGate::new(),
        ),
        Err(ChtError::InvalidInput {
            field: "unsteady.implicit_buoyancy",
            ..
        })
    ));
}
