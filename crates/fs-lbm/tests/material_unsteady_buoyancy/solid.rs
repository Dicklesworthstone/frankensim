//! Independent caloric/finite-volume oracles for sourced solid transients.

use super::*;
use fs_evidence::ValidityDomain;
use fs_lbm::material::conjugate::buoyant::SolidThermalConfig;
use fs_lbm::material::conjugate::solid::{SolidConductivityLaw, SolidThermalLaw};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point,
};
use fs_material::thermal::transport::storage::HeatStorageAmount;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec, Temperature};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn card(capacity_slope: f64, conductivity_slope: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        // Evolving source density MUST NOT change the fixed reference mass.
        ("rho", Density::DIMS, [3.0, 3.0, 2.97, 2.7]),
        (
            "cp",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            [
                2.0,
                2.0,
                2.0 + 10.0 * capacity_slope,
                2.0 + 190.0 * capacity_slope,
            ],
        ),
        (
            "k",
            THERMAL_CONDUCTIVITY_DIMS,
            [
                1.0,
                1.0,
                1.0 + 10.0 * conductivity_slope,
                1.0 + 100.0 * conductivity_slope,
            ],
        ),
    ] {
        let response = if name == "k" {
            quantity(QuantityKind::ThermalConductivity)
        } else {
            QuantitySpec::dimensional(dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, response),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: [290.0, 300.0, 310.0, 400.0]
                        .into_iter()
                        .zip(values)
                        .collect(),
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        quantity(QuantityKind::AbsoluteTemperature),
                        290.0,
                        400.0,
                    )
                    .with_quantity(
                        "pressure",
                        quantity(QuantityKind::Pressure),
                        100000.0,
                        100000.0,
                    ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic transient solid calorimetry, not measurements".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic solid".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn law(card: &MaterialCard, span: [f64; 2]) -> SolidThermalLaw<'_> {
    let point = QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 300.0)
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100000.0)
        .unwrap();
    let requirements = [
        ("rho", QuantitySpec::dimensional(Density::DIMS)),
        ("cp", QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS)),
        ("k", quantity(QuantityKind::ThermalConductivity)),
    ]
    .map(|(name, response)| {
        ScalarPropertyRequirement::try_with_quantity(
            name,
            response,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()
    });
    let source = resolve_material_state_point(
        card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    SolidConductivityLaw::try_new(card, source, "T", ["k"; 3])
        .unwrap()
        .with_heat_storage("rho", "cp", span)
        .unwrap()
}

// Closed-form integral of the SYNTHETIC piecewise-affine Cp source, from 300 K.
// Keeping this independent of PreparedHeatStorage catches missed source knots.
fn specific_energy(temperature: f64, slope: f64) -> f64 {
    let theta = temperature - 300.0;
    if theta <= 0.0 {
        2.0 * theta
    } else if theta <= 10.0 {
        2.0 * theta + 0.5 * slope * theta * theta
    } else {
        let beyond = theta - 10.0;
        20.0 + 50.0 * slope + (2.0 + 10.0 * slope) * beyond + slope * beyond * beyond
    }
}

fn pair() -> VoxelDomain {
    let mut domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    domain.set(0, 0, 0, Voxel::Solid(0));
    domain
}

#[test]
fn g1_exact_piecewise_cp_and_endpoint_k_close_native_two_cell_heat_balance() {
    let domain = pair();
    let reference = FluidTransportState::try_new(1.0, 20.0, 2.0, 1.0).unwrap();
    let fluid = MaterialChtFluid::try_new(&reference).unwrap();
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let unsteady = UnsteadyConfig::new(2.0, 3);
    let mut temperatures = Vec::new();
    for (cp_slope, k_slope) in [(0.0, 0.0), (0.2, 0.0), (0.2, 0.01)] {
        let card = card(cp_slope, k_slope);
        let solid = law(&card, [290.0, 400.0]);
        let result = fluid
            .march_natural_convection_with_material_solids(
                &domain,
                std::slice::from_ref(&solid),
                &flow,
                &thermal,
                &EnergyConfig::default(),
                &[300.0, 400.0],
                limits([0.0; 3]),
                &unsteady,
                SolidThermalConfig::default(),
                |_| 1.0,
                |_| 1.0,
                None,
                &CancelGate::new(),
                |_, _, _| Ok(reference),
            )
            .unwrap();
        let (mut ts, mut tf) = (300.0, 400.0);
        for (index, record) in result.fluid.march.energy_records.iter().enumerate() {
            // Harmonic two half-cell resistance; the implicit fluid equation
            // eliminates Tf and leaves this independent scalar solid balance.
            let (mut low, mut high) = (ts, tf);
            for _ in 0..100 {
                let mid = 0.5 * (low + high);
                let k = 1.0 + k_slope * (mid - 300.0);
                let conductance = 1.0 / (0.5 / k + 0.5 / 2.0);
                let heat = 2.0 * conductance * 20.0 * (tf - mid) / (20.0 + 2.0 * conductance);
                let stored = 3.0 * (specific_energy(mid, cp_slope) - specific_energy(ts, cp_slope));
                if stored > heat {
                    high = mid;
                } else {
                    low = mid;
                }
            }
            let next_s = 0.5 * (low + high);
            let stored = 3.0 * (specific_energy(next_s, cp_slope) - specific_energy(ts, cp_slope));
            tf -= stored / 20.0;
            ts = next_s;
            close(
                result.solid_records[index].solid_stored_energy_change_j,
                stored,
                2e-9,
            );
            assert!(record.stored_energy_change_j.abs() < 2e-7, "{record:?}");
            assert!(record.closure_j.abs() < 2e-7, "{record:?}");
            assert_eq!(record.source_j, 0.0);
            assert_eq!(record.boundary_outflow_j, 0.0);
        }
        close(result.fluid.march.temperature[0], ts, 1e-10);
        close(result.fluid.march.temperature[1], tf, 1e-10);
        assert!(ts > 310.0); // The actual trajectory crosses the Cp slope break.
        close(
            result.solid_energy_j[0],
            3.0 * (20.0 + specific_energy(ts, cp_slope)),
            1e-10,
        );
        assert_eq!(result.solid_energy_j[1], 0.0);
        assert_eq!(
            solid.storage().components()[0].amount(),
            HeatStorageAmount::MassKg(3.0)
        );
        let endpoint = &result.solid_endpoint[0];
        assert_eq!(endpoint.temperature_k, result.fluid.march.temperature[0]);
        assert_eq!(endpoint.material.card_identity(), card.content_hash());
        assert!(endpoint.material.property("rho").unwrap().value_si() < 3.0);
        for name in ["rho", "cp", "k"] {
            let receipt = &endpoint.material.property(name).unwrap().answer().receipt;
            assert_eq!(
                receipt.selected,
                solid
                    .conductivity()
                    .reference()
                    .property(name)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
            assert!(receipt.query_point.contains(&("pressure".into(), 100000.0)));
        }
        close(
            endpoint.conductivity_w_m_k[0],
            1.0 + k_slope * (ts - 300.0),
            1e-12,
        );
        temperatures.push(ts);
    }
    assert!(temperatures[1] < temperatures[0]); // More heat capacity delays warming.
    assert!(temperatures[2] > temperatures[1]); // Higher actual k raises heat transfer.
}

#[test]
fn g1_fixed_temperature_storage_and_radiation_use_true_source_integral_once() {
    let dx = 0.01;
    let domain = VoxelDomain::from_fn(1, 1, 2, dx, |p| {
        if p[2] < dx {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let reference = FluidTransportState::try_new(1.0, 1.0, 1e-12, 1.0).unwrap();
    let fluid = MaterialChtFluid::try_new(&reference).unwrap();
    let card = card(0.2, 0.01);
    let solid = law(&card, [290.0, 400.0]);
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.fixed_temperature = vec![(0, 330.0), (0, 350.0)];
    let mut surroundings = [Some(300.0); 6];
    surroundings[4] = None;
    let radiation = RadiationConfig::new(vec![0.9], surroundings);
    let result = fluid
        .march_natural_convection_with_material_solids(
            &domain,
            &[solid],
            &SimpleConfig::new([FvBoundary::wall(); 6]),
            &thermal,
            &EnergyConfig::default(),
            &[300.0; 2],
            limits([0.0; 3]),
            &UnsteadyConfig::new(2.0, 2),
            SolidThermalConfig::default(),
            |_| 1.0,
            |_| 1.0,
            Some(&radiation),
            &CancelGate::new(),
            |_, _, _| Ok(reference),
        )
        .unwrap();
    let expected = 3.0 * domain.volume(0) * specific_energy(350.0, 0.2);
    close(
        result.solid_records[0].solid_stored_energy_change_j,
        expected,
        1e-13,
    );
    close(
        result.solid_records[1].solid_stored_energy_change_j,
        0.0,
        1e-13,
    );
    let radiated_w = 0.9 * STEFAN_BOLTZMANN * dx * dx * (350.0_f64.powi(4) - 300.0_f64.powi(4));
    for record in &result.fluid.march.energy_records {
        close(record.radiated_w, radiated_w, 1e-11);
        close(record.closure_j, 0.0, 1e-10);
        close(
            record.boundary_outflow_j,
            -record.stored_energy_change_j,
            1e-10,
        );
    }
    close(
        result.fluid.march.energy_records[0].stored_energy_change_j,
        expected,
        1e-10,
    );
}

#[test]
fn g1_large_predictor_backtracks_without_clamping_the_accepted_source_state() {
    let domain = pair();
    let reference = FluidTransportState::try_new(1.0, 2.0, 1e-12, 1.0).unwrap();
    let fluid = MaterialChtFluid::try_new(&reference).unwrap();
    let card = card(2.0, 0.0);
    let solid = law(&card, [290.0, 310.0]);
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.power_w = vec![90.0, 0.0];
    thermal.fixed_temperature = vec![(1, 300.0)];
    let result = fluid
        .march_natural_convection_with_material_solids(
            &domain,
            &[solid],
            &SimpleConfig::new([FvBoundary::wall(); 6]),
            &thermal,
            &EnergyConfig::default(),
            &[300.0; 2],
            limits([0.0; 3]),
            &UnsteadyConfig::new(1.0, 1),
            SolidThermalConfig::default(),
            |_| 1.0,
            |_| 1.0,
            None,
            &CancelGate::new(),
            |_, _, _| Ok(reference),
        )
        .unwrap();
    // The frozen-reference-capacity predictor is 315 K, outside explicit Cp
    // support. The true integral gives 3*(2*theta+theta^2)=90 J, inside it.
    let expected = 299.0 + 31.0_f64.sqrt();
    close(result.fluid.march.temperature[0], expected, 1e-10);
    assert!(result.solid_records[0].iterations > 2);
    close(
        result.solid_records[0].solid_stored_energy_change_j,
        90.0,
        2e-10,
    );
    close(result.fluid.march.energy_records[0].closure_j, 0.0, 2e-8);
    assert!(result.solid_endpoint[0].temperature_k < 310.0);
}

#[test]
fn g0_source_span_initial_endpoint_budget_and_cancellation_refuse_atomically() {
    let domain = pair();
    let reference = FluidTransportState::try_new(1.0, 2.0, 1e-12, 1.0).unwrap();
    let fluid = MaterialChtFluid::try_new(&reference).unwrap();
    let card = card(2.0, 0.0);
    let solid = law(&card, [290.0, 310.0]);
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.power_w = vec![90.0, 0.0];
    thermal.fixed_temperature = vec![(1, 300.0)];
    let unsteady = UnsteadyConfig::new(1.0, 1);
    let run = |setup: &ThermalSetup, initial: &[f64], controls, gate: &CancelGate| {
        fluid.march_natural_convection_with_material_solids(
            &domain,
            std::slice::from_ref(&solid),
            &flow,
            setup,
            &EnergyConfig::default(),
            initial,
            limits([0.0; 3]),
            &unsteady,
            controls,
            |_| 1.0,
            |_| 1.0,
            None,
            gate,
            |_, _, _| Ok(reference),
        )
    };
    assert!(matches!(
        run(
            &thermal,
            &[311.0, 300.0],
            SolidThermalConfig::default(),
            &CancelGate::new()
        ),
        Err(MaterialBuoyancyError::SolidThermal { cell: 0, .. })
    ));
    let mut unsupported = thermal.clone();
    unsupported.fixed_temperature.push((0, 311.0));
    assert!(matches!(
        run(
            &unsupported,
            &[300.0; 2],
            SolidThermalConfig::default(),
            &CancelGate::new()
        ),
        Err(MaterialBuoyancyError::SolidThermal { cell: 0, .. })
    ));
    let mut overheat = thermal.clone();
    overheat.power_w[0] = 400.0;
    assert!(
        run(
            &overheat,
            &[300.0; 2],
            SolidThermalConfig::default(),
            &CancelGate::new()
        )
        .is_err()
    );
    let mut low_budget = SolidThermalConfig::default();
    low_budget.maximum_iterations = 1;
    assert!(run(&thermal, &[300.0; 2], low_budget, &CancelGate::new()).is_err());
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        run(&thermal, &[300.0; 2], SolidThermalConfig::default(), &gate),
        Err(MaterialBuoyancyError::Solver(ChtError::Cancelled))
    ));
    // A refusal changed no caller-owned source or setup: the same input retries.
    assert!(
        run(
            &thermal,
            &[300.0; 2],
            SolidThermalConfig::default(),
            &CancelGate::new()
        )
        .is_ok()
    );
    assert_eq!(solid.storage().temperature_span_k(), [290.0, 310.0]);
    assert_eq!(solid.conductivity().reference().query().axes()["T"], 300.0);
    assert_eq!(thermal.power_w, vec![90.0, 0.0]);
}
