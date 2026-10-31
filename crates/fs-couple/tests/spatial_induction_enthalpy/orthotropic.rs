//! G1/G3: oriented electrical sources share the native 3D field/enthalpy residual.
//! Synthetic curves qualify equations, not measured material data.

use super::*;
use fs_matdb::{InterpolationPolicy, PropertyValue};
use fs_material::{
    conductor::spatial::{ResolvedOrthotropicConductor3, SpatialConductorMaterial},
    state_point::{MaterialPropertySelection, resolve_material_state_point},
};
use fs_qty::Temperature;

fn keys() -> [PropertyKey; 3] {
    ["sigma-longitudinal", "sigma-transverse", "sigma-through"]
        .map(|name| PropertyKey::new(name, ElectricalPropertyBasis::Conductivity.dimensions()))
}

fn card(beta: f64) -> MaterialCard {
    let base = source::absolute_affine_magnetic_card(
        1.0, 0.0, 1.0, 0.0, ElectricalPropertyBasis::Conductivity,
    );
    let mut claims = base.claims().clone();
    let template = base.claims().claims_for("electrical_conductivity")[0].1.clone();
    for (key, sigma) in keys().into_iter().zip([2.0, 5.0, 11.0]) {
        let mut claim = template.clone();
        claim.key = key;
        claim.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: [200.0, 300.0, 2000.0]
                .map(|t| (t, sigma * (1.0 + beta * (t - 300.0))))
                .to_vec(),
            dims: ElectricalPropertyBasis::Conductivity.dimensions(),
        };
        claim.interpolation = InterpolationPolicy::LinearInside;
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

fn problem(card: &MaterialCard) -> SpatialVectorMagnetostatics {
    let basis = ElectricalPropertyBasis::Conductivity;
    let original = super::problem(card, basis);
    let context = source::absolute_point(300.0);
    let scalar = source::bundle_at(card, &context, basis);
    let mut requirements = scalar.properties().iter()
        .map(|property| property.requirement().clone()).collect::<Vec<_>>();
    requirements.extend(keys().iter().map(|key| basis.requirement_for(key).unwrap()));
    let material = resolve_material_state_point(
        card, &context, &requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    SpatialVectorMagnetostatics::new(
        original.positions_m().to_vec(),
        original.mesh().tets.clone(),
        vec![SpatialMagneticMaterial::Linear(material); 6],
        assembly(),
        &mut || false,
    ).unwrap()
}

fn initial(
    problem: &SpatialVectorMagnetostatics,
    axes: [[f64; 3]; 3],
    amplitude: f64,
) -> SpatialInduction<'_, ResolvedOrthotropicConductor3> {
    let boundary = VectorMagneticBoundary {
        potential: problem.exterior_faces().into_iter().map(|face| VectorPotentialFace {
            face,
            edge_integrals_wb: [0.0; 3],
        }).collect(),
    };
    let solve = CheckedCgConfig {
        tolerance: 1e-12,
        max_iterations: 1000,
        max_corrections: 2,
    };
    let initial = problem.solve(
        &vec![[8.0 * amplitude; 3]; 6],
        &boundary,
        VectorMagneticControls {
            field: solve,
            gauge: solve,
            maximum_gauge_product_terms: 100000,
            current_tolerance_a: 1e-9,
            flux_tolerance_wb: 1e-10,
            energy_tolerance_j: 1e-9,
        },
        &mut || false,
    ).unwrap();
    let conductors = problem.materials().iter().map(|material| {
        let SpatialMagneticMaterial::Linear(material) = material else { unreachable!() };
        Some(ResolvedOrthotropicConductor3::try_from_material(
            material.clone(), ElectricalPropertyBasis::Conductivity, keys(), axes,
        ).unwrap())
    }).collect();
    SpatialInduction::new_with_conductors(
        &initial,
        conductors,
        SpatialInductionLimits {
            maximum_unknowns: 1000,
            maximum_conducting_cells: 100,
            maximum_assembly_terms: 10000,
        },
        1e-10,
        &mut || false,
    ).unwrap()
}

#[test]
fn rotated_conductivity_and_temperature_close_the_same_field_heat_history() {
    let c = core::f64::consts::FRAC_1_SQRT_2;
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let positive = [[c, c, 0.0], [-c, c, 0.0], [0.0, 0.0, 1.0]];
    let negative = [[c, -c, 0.0], [c, c, 0.0], [0.0, 0.0, 1.0]];
    let dt = 0.1;
    let a0 = 3.0;
    let mut endpoints = Vec::new();
    for (axes, beta, m0) in [
        (identity, 0.0, 0.45),
        (identity, 0.005, 0.45),
        (positive, 0.005, 0.4625),
        (negative, 0.005, 0.4375),
    ] {
        let card = card(beta);
        let curve = curve(&card);
        let thermal = thermal(&curve, 0.0);
        let problem = problem(&card);
        let field = initial(&problem, axes, a0);
        let context = source::absolute_point(300.0);
        let cells = [MaterialCell {
            cell: CellId(7),
            material_id: MaterialId(9),
            curve: curve.curve(),
            point: ThermalMaterialPoint {
                material: &card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        }];
        let regions = regions();
        let bindings = laws(&card, &problem, 6.1);
        let network = SpatialMagneticThermalNetwork {
            cells: &cells,
            links: &[],
            magnetic_regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        let candidate = network.step_candidate(
            &thermal, dt, &field, &[[0.0; 3]; 6], &[[0.0; 3]; 8],
            &bindings, thermal_controls(), field_controls(), assembly(), &[0.0], || false,
        ).unwrap();

        // Independent exact elimination on the six-tet cube, with mu=1:
        // K=4; M=tr(sigma)/40-(sigma_xy+sigma_xz+sigma_yz)/120.
        // Principal sources share an authored affine temperature scale.
        // The accepted old impressed current ramps to zero over this step.
        let mass = |t: f64| m0 * (1.0 + beta * (t - 300.0));
        let amplitude = |t: f64| a0 * mass(t) / (mass(t) + 2.0 * dt);
        let heat = |t: f64| mass(t) * (amplitude(t) - a0).powi(2) / dt;
        let expected = bracket(300.0, 310.0, |t| 4.0 * (t - 300.0) - heat(t));
        let temperature = candidate.thermal().temperature_k(CellId(7)).unwrap();
        let actual_amplitude = candidate.field().potential_wb()[diagonal(&problem)];
        near(temperature, expected, 2e-8);
        near(actual_amplitude, amplitude(expected), 2e-9);
        near(candidate.report().coupled.electrical.heat_j, heat(expected), 2e-8);
        near(candidate.report().coupled.deposited_joule_j, heat(expected), 2e-8);
        near(candidate.report().coupled.thermal.stored_energy_change_j, heat(expected), 2e-8);
        near(candidate.report().external_heat_j, 0.0, 2e-8);
        near(candidate.report().total_energy_residual_j, 0.0, 2e-8);
        for source in candidate.field().conductors().iter().flatten() {
            assert_eq!(source.principal_axes(), Some(axes));
            for (key, sigma) in keys().iter().zip([2.0, 5.0, 11.0]) {
                near(
                    source.material().property_by_key(key).unwrap().value_si(),
                    sigma * (1.0 + beta * (temperature - 300.0)),
                    2e-10,
                );
            }
        }
        assert_eq!(field.time_s(), 0.0);
        near(thermal.temperature_k(CellId(7)).unwrap(), 300.0, 1e-12);
        let continued = network.step_candidate(
            candidate.thermal(), 2.0 * dt, candidate.field(),
            &[[0.0; 3]; 6], &[[0.0; 3]; 8], &bindings,
            thermal_controls(), field_controls(), assembly(), &[0.0], || false,
        ).unwrap();
        near(continued.report().total_energy_residual_j, 0.0, 2e-8);
        assert_eq!(continued.field().time_s(), 2.0 * dt);
        assert!(matches!(
            network.step_candidate(
                &thermal, dt, &field, &[[0.0; 3]; 6], &[[0.0; 3]; 8],
                &bindings, thermal_controls(), field_controls(), assembly(), &[0.0], || true,
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        endpoints.push(temperature);
    }
    assert!((endpoints[0] - endpoints[1]).abs() > 1e-6);
    assert!((endpoints[2] - endpoints[3]).abs() > 1e-6);
}
