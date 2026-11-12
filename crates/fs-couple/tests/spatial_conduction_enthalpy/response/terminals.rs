//! Actual-load electrode closure stays inside the original nonlinear source domain.
use super::*;
use fs_em::conduction::response::terminals::{
    ConductionTerminalControl, ConductionTerminalDrive, ConductionTerminalLimits,
};

#[test]
fn narrow_current_domain_and_floating_electrode_close_heat_and_continue_under_huge_bias() {
    let coordinate = ConductorMagnitude::CurrentDensity;
    let r = [1.0, 2.0];
    let beta = [0.2, 0.08];
    let cards = [0, 1].map(|i| {
        let (card, _) = nonlinear_material(r[i], CP[i], ALPHA[i], beta[i], coordinate);
        let mut claims = ClaimSet::new();
        for (_, original) in card.claims().claims_ordered() {
            let mut claim = original.clone();
            if claim.key.name() == "original_nonlinear_transport" {
                // The authored operating range excludes artificial unit-voltage
                // probes but includes every actual J=2 A/m² endpoint.
                claim.validity =
                    claim
                        .validity
                        .with_quantity("operating", coordinate.quantity(), 1.5, 2.5);
                // A claim's validity must enclose every source sample of its surface (fs-matdb
                // insert_claim), so the surface is cut to the same range. It is affine in T and
                // the current density, so every value inside [1.5, 2.5] is unchanged.
                let basis = ElectricalPropertyBasis::Resistivity;
                claim.value = PropertyValue::Surface2D {
                    axes: ["T".into(), "operating".into()],
                    axis_dims: [Temperature::DIMS, coordinate.quantity().dims()],
                    samples: [[200.0, 1.5], [2000.0, 1.5], [200.0, 2.5], [2000.0, 2.5]]
                        .into_iter()
                        .map(|p| (p, r[i] + ALPHA[i] * (p[0] - 300.0) + beta[i] * p[1]))
                        .collect(),
                    triangles: vec![[0, 1, 2], [1, 3, 2]],
                    dims: basis.dimensions(),
                };
            }
            claims.insert_claim(claim).unwrap();
        }
        MaterialCard::assemble(card.id().clone(), claims, vec![]).unwrap()
    });
    let context_at = |t| {
        point(t)
            .with_quantity("operating", coordinate.quantity(), 1.75)
            .unwrap()
    };
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let curves = cards.each_ref().map(|card| {
        ResolvedHeatCapacityCurve::try_new(
            card,
            HeatCapacitySourceQuery {
                heat_capacity: &heat,
                density: &density,
                temperature_axis: "T",
                lower: &context_at(200.0),
                upper: &context_at(2000.0),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            0.0,
            1e-6,
        )
        .unwrap()
    });
    let context = context_at(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &key);
    let regions = regions();
    let layout = SpatialConductionThermalNetwork {
        cells: &cells,
        links: &LINKS,
        regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let original = problem_at(&cards, coordinate, &context);
    let (mut positions, mut tets) = slab(&[0.0, 1.0]);
    let (second, second_tets) = slab(&[1.0, 2.0]);
    positions.extend(second);
    tets.extend(second_tets.into_iter().map(|tet| tet.map(|v| v + 8)));
    let mut problem = SpatialConductionResponse::new(
        positions,
        tets,
        original.sources().to_vec(),
        limits(),
        &mut || false,
    )
    .unwrap()
    .with_relaxation(0.5)
    .unwrap();
    let electrodes = [0.0, 1.0, 2.0].map(|x| {
        problem
            .exterior_faces()
            .filter(|&face| {
                problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0]) == [x; 3]
            })
            .collect::<Vec<_>>()
    });
    let terminals = electrodes.each_ref().map(Vec::as_slice);
    let drive = ConductionTerminalDrive {
        terminals: &terminals,
        control: ConductionTerminalControl::Mixed {
            current_terminals: &[2, 1],
            currents_a: &[2.0, 0.0],
            voltages_v: &[1e250],
        },
        limits: ConductionTerminalLimits {
            maximum_terminals: 3,
            maximum_workspace_values: 20_000,
        },
    };
    let mut thermal = initial(&curves, 999.0, None);
    let mut old = [0.0; 2];
    let mut solve = config();
    solve.maximum_circuit_trials = 4096;
    for target_s in [0.1, 0.2] {
        let mut polls = 0usize;
        let next = layout
            .step_response_candidate(
                &thermal,
                target_s,
                &problem,
                drive,
                limits(),
                field_controls(),
                controls(),
                solve,
                &[0.0; 2],
                || {
                    polls += 1;
                    false
                },
            )
            .unwrap();
        let (theta, e) = endpoint(2.0, 0.1, old, r, beta, false);
        let step = &next.report().electrical;
        let observed = step.field.terminals.as_ref().unwrap();
        for (i, cell) in [CellId(30), CellId(10)].into_iter().enumerate() {
            near(
                next.thermal().temperature_k(cell).unwrap(),
                300.0 + theta[i],
            );
        }
        near(observed.voltage_v[0], 0.0);
        near(observed.voltage_v[1], e[0]);
        near(observed.voltage_v[2], e[0] + e[1]);
        for (actual, expected) in observed.current_a.iter().zip([-2.0, 0.0, 2.0]) {
            near(*actual, expected);
        }
        // Absolute observations round to the huge common bias; the retained
        // report must still carry the finite physical voltage and source work.
        assert!(step.field.field.voltage_v.iter().all(|&v| v == 1e250));
        near(step.supplied_work_j, 0.2 * (e[0] + e[1]));
        near(step.heat_j, step.supplied_work_j);
        near(next.thermal().energy_j() - thermal.energy_j(), step.heat_j);
        near(
            CP[0] * (theta[0] - old[0]) + CP[1] * (theta[1] - old[1]),
            step.heat_j,
        );
        near(next.report().total_energy_residual_j, 0.0);
        assert_eq!(step.field.field_solves, step.field.iterations);
        for source in next.problem().conductors() {
            near(source.material().query().axes()["operating"], 2.0);
        }
        let mut bounded = solve;
        bounded.maximum_circuit_trials = next.report().circuit_trials - 1;
        assert!(matches!(
            layout.step_response_candidate(
                &thermal,
                target_s,
                &problem,
                drive,
                limits(),
                field_controls(),
                controls(),
                bounded,
                &[0.0; 2],
                || false,
            ),
            Err(CircuitEnthalpyError::CircuitBudget)
        ));
        let mut stop = 0usize;
        assert!(matches!(
            layout.step_response_candidate(
                &thermal,
                target_s,
                &problem,
                drive,
                limits(),
                field_controls(),
                controls(),
                solve,
                &[0.0; 2],
                || {
                    stop += 1;
                    stop >= polls - 1
                },
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(thermal.time_s(), target_s - 0.1);
        (thermal, problem, _) = next.into_parts();
        old = theta;
    }
    assert_eq!(context.axes()["operating"], 1.75);
    assert_eq!(
        original.conductors()[0].material().query().axes()["operating"],
        1.75
    );
}
