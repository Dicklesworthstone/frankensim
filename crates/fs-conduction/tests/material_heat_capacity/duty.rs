//! Material heat storage reaches scheduled-power diffusion and source support.

use super::*;
use fs_conduction::duty::{
    DutyCycle, DutySegment, march_duty_cycle_with_capacity,
};
use fs_conduction::transient::{HeatCapacity, TransientConfig, VolumetricHeatCapacity};

#[test]
fn g1_two_material_ramp_heating_retains_each_caloric_mass_and_diffusion_rate() {
    let mesh = mesh(true);
    let first = card(2.0, 10.0, 2.0, false);
    let second = card(3.0, 20.0, 4.0, false);
    let ca = capacity(&first, 360.0);
    let cb = capacity(&second, 360.0);
    let ka = ca
        .isotropic_conductivity(&first, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    let kb = cb
        .isotropic_conductivity(&second, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    let assigned = ElementMaterials::new(
        MaterialTable::new([(MaterialId(7), ka), (MaterialId(2), kb)]).unwrap(),
        vec![MaterialId(7), MaterialId(2)],
    )
    .unwrap();
    // Deliberately reverse source-table order relative to element order.
    let capacities = [(MaterialId(2), cb), (MaterialId(7), ca)];
    let fallback = ConductivityModel::isotropic_declared(99.0).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let source = ScalarField::Uniform(30.0);
    let initial = [300.0, 310.0, 290.0, 300.0, 300.0, 310.0, 290.0, 300.0];
    let cycle = DutyCycle::new(vec![DutySegment::ramp(0.75, 0.0, 1.0).unwrap()]).unwrap();

    for theta in [1.0, 0.75, 0.5] {
        let config = TransientConfig::new(theta, 0.25, controls().linear).unwrap();
        let solved = with_cx(|cx| {
            march_duty_cycle_with_capacity(
                cx,
                ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &fallback,
                    source: &source,
                    element_materials: Some(&assigned),
                },
                HeatCapacity::MaterialById(&capacities),
                &config,
                &cycle,
                10.0, // two tetrahedra, each with volume 1/6 m^3
                &initial,
                f64::INFINITY,
            )
            .unwrap()
        });

        // Independent P1 mode [0,1,-1,0] on the unit right tetrahedron:
        // eigenvalue 4*k/(rho*Cp). Uniform heating affects only the mean.
        let scale_seconds = 0.375 + (theta - 0.5) * 0.25;
        let mut stored_j = 0.0;
        for (element, (volumetric, conductivity)) in [(20.0, 2.0), (60.0, 4.0)]
            .into_iter()
            .enumerate()
        {
            let z = 4.0 * conductivity * 0.25 / volumetric;
            let ratio = (1.0 - (1.0 - theta) * z) / (1.0 + theta * z);
            let amplitude = 10.0 * ratio.powi(3);
            let mean = 300.0 + 30.0 * scale_seconds / volumetric;
            let expected = [mean, mean + amplitude, mean - amplitude, mean];
            for (local, &temperature) in expected.iter().enumerate() {
                let vertex = 4 * element + local;
                assert!(
                    (solved.temperature[vertex] - temperature).abs() < 1e-8,
                    "theta={theta} material={element} vertex={local}: {} != {temperature}",
                    solved.temperature[vertex]
                );
                stored_j += volumetric / 24.0 * (solved.temperature[vertex] - initial[vertex]);
            }
        }
        let supplied_j = 10.0 * scale_seconds;
        assert!((stored_j - supplied_j).abs() < 1e-8);
        assert!((solved.energy.delivered_j() - supplied_j).abs() < 1e-12);
    }
}

#[test]
fn g1_inter_material_diffusion_conserves_heat_and_reaches_the_capacity_weighted_equilibrium() {
    let points = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ];
    let mesh = ConductionMesh::new(
        TetComplex::from_tets(points.len(), vec![[0, 1, 2, 3], [0, 2, 1, 4]]),
        points,
    )
    .unwrap();
    let first = ConductivityModel::isotropic_declared(2.0).unwrap();
    let second = ConductivityModel::isotropic_declared(4.0).unwrap();
    let fallback = first.clone();
    let assigned = ElementMaterials::new(
        MaterialTable::new([(MaterialId(7), first), (MaterialId(2), second)]).unwrap(),
        vec![MaterialId(7), MaterialId(2)],
    )
    .unwrap();
    let capacities = [
        VolumetricHeatCapacity::declared(20.0).unwrap(),
        VolumetricHeatCapacity::declared(60.0).unwrap(),
    ];
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let source = ScalarField::Uniform(0.0);
    let initial = [300.0, 300.0, 300.0, 340.0, 280.0];
    let cycle = DutyCycle::new(vec![DutySegment::constant(200.0, 1.0).unwrap()]).unwrap();
    let config = TransientConfig::backward_euler(2.0, controls().linear).unwrap();
    let solved = with_cx(|cx| {
        march_duty_cycle_with_capacity(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &fallback,
                source: &source,
                element_materials: Some(&assigned),
            },
            HeatCapacity::PerElement(&capacities),
            &config,
            &cycle,
            0.0,
            &initial,
            f64::INFINITY,
        )
        .unwrap()
    });
    // Each tet contributes rho*Cp*V/4 to its incident vertices. The common
    // face stores both materials' contributions, rather than an averaged Cp.
    let weights = [80.0 / 24.0, 80.0 / 24.0, 80.0 / 24.0, 20.0 / 24.0, 60.0 / 24.0];
    let equilibrium = 298.75;
    let mut stored_j = 0.0;
    for (vertex, &temperature) in solved.temperature.iter().enumerate() {
        assert!((temperature - equilibrium).abs() < 1e-8);
        stored_j += weights[vertex] * (temperature - initial[vertex]);
    }
    assert!(stored_j.abs() < 1e-8, "internal diffusion created {stored_j} J");
    assert_eq!(solved.energy.delivered_j(), 0.0);
}

#[test]
fn g4_duty_source_domains_cover_initial_prescribed_and_all_accepted_states() {
    let mesh = mesh(false);
    let original = card(2.0, 10.0, 2.0, false);
    let narrow = capacity(&original, 300.5);
    let wide = capacity(&original, 360.0);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let narrow_k = ConductivityModel::isotropic(
        fs_conduction::ConductivityTable::declared_curve(vec![(280.0, 2.0), (300.2, 2.0)])
            .unwrap(),
    );
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let prescribed = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "unsupported wall",
            |face| face.centroid[0] == 0.0,
            ThermalBc::dirichlet(301.0).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let source = ScalarField::Uniform(30.0);
    let initial = [300.0; 4];
    let short = DutyCycle::new(vec![DutySegment::constant(0.25, 1.0).unwrap()]).unwrap();
    let long = DutyCycle::new(vec![DutySegment::constant(0.5, 1.0).unwrap()]).unwrap();
    let config = TransientConfig::crank_nicolson(0.25, controls().linear).unwrap();
    with_cx(|cx| {
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            source: &source,
            element_materials: None,
        };
        let accepted = march_duty_cycle_with_capacity(
            cx, problem, HeatCapacity::MaterialUniform(&narrow),
            &config, &short, 5.0, &initial, f64::INFINITY,
        )
        .unwrap();
        assert!(accepted.temperature.iter().all(|t| (t - 300.375).abs() < 1e-9));
        for (candidate, storage, cycle, start) in [
            (problem, &narrow, &long, initial),
            (ConductionProblem { material: &narrow_k, ..problem }, &wide, &short, initial),
            (ConductionProblem { boundary: &prescribed, ..problem }, &narrow, &short, initial),
            (problem, &narrow, &short, [300.6; 4]),
        ] {
            assert!(matches!(
                march_duty_cycle_with_capacity(
                    cx, candidate, HeatCapacity::MaterialUniform(storage),
                    &config, cycle, 5.0, &start, f64::INFINITY,
                ),
                Err(ConductionError::OutsideTemperatureSpan { .. })
            ));
        }
        assert_eq!(initial, [300.0; 4]);
    });
}
