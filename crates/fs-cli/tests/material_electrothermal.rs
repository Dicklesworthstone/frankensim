//! Actual material-card -> constitutive owner -> spatial electrical/thermal
//! integration. The dimensioned data below are synthetic analytic fixtures,
//! NOT measurements of a real metal or a physical validation claim.

use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_lbm::conjugate::material_coupling::{
    MaterialCoupledSolution, MaterialCouplingConfig, solve_material_electrothermal,
};
use fs_lbm::conjugate::{
    ChtError, Conductor, ElectricSetup, Electrode, ElectrodeDrive, EnergyConfig, FlowField,
    FluidProperties, SolidMaterial, ThermalFace, ThermalSetup, Voxel, VoxelDomain, solve_energy,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::conductor::{
    ConductorError, ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY,
    resolve_conductor_state, resolve_uniform_conductor,
};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
    resolve_material_state_point,
};
use fs_qty::Dims;

fn card(name: &str, resistivity_scale: f64, conductivity: f64, maximum_k: f64) -> MaterialCard {
    let provenance = Provenance {
        source: "authored analytic electrothermal fixture; no observation claim".into(),
        license: "same as repository (synthetic test data)".into(),
        artifact: None,
    };
    let validity = ValidityDomain::unconstrained()
        .with("T", 300.0, maximum_k)
        .with("pressure", 100_000.0, 110_000.0);
    let uncertainty = UncertaintyModel::RelativeHalfWidth {
        fraction: 0.01,
        confidence: 0.95,
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Dims([0, 0, 0, 1, 0, 0]),
                knots: vec![
                    (300.0, 2.0 * resistivity_scale),
                    (
                        maximum_k,
                        2.0 * resistivity_scale * (1.0 + 0.05 * (maximum_k - 300.0)),
                    ),
                ],
                dims: ELECTRICAL_RESISTIVITY_DIMS,
            },
            validity: validity.clone(),
            uncertainty: uncertainty.clone(),
            interpolation: InterpolationPolicy::LinearInside,
            observations: Vec::new(),
            provenance: provenance.clone(),
        })
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: conductivity,
                dims: THERMAL_CONDUCTIVITY_DIMS,
            },
            validity,
            uncertainty,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance,
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "synthetic homogeneous solid".into(),
            process: "analytic fixture".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

fn point(temperature: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", temperature)
        .unwrap()
        .with("pressure", pressure)
        .unwrap()
}

fn source_error(field: &'static str, error: impl core::fmt::Display) -> ChtError {
    ChtError::InvalidInput {
        field,
        reason: error.to_string(),
    }
}

struct Run {
    result: MaterialCoupledSolution,
    electrical_states: Vec<ResolvedMaterialStatePoint>,
    thermal_states: Vec<ResolvedMaterialStatePoint>,
}

/// No branch on chemistry, material name or selected variant. Electrical
/// coefficients come from source curves; the thermal fixtures declare constant
/// k within their support. Numerical kernels and geometry stay unchanged.
fn run(card: &MaterialCard) -> Result<Run, ChtError> {
    let gate = CancelGate::new();
    let domain = VoxelDomain::from_fn(3, 1, 1, 1.0, |_| Voxel::Solid(0))?;
    let initial = resolve_conductor_state(
        card,
        &point(300.0, 101_325.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .map_err(|e| source_error("material.electrical_resistivity", e))?;
    let electric = ElectricSetup {
        conductors: vec![Some(Conductor::new(initial.resistivity_ohm_m()))],
        electrodes: vec![
            Electrode {
                lo: [0, 0, 0],
                hi: [1, 1, 1],
                drive: ElectrodeDrive::Voltage(0.0),
            },
            Electrode {
                lo: [2, 0, 0],
                hi: [3, 1, 1],
                drive: ElectrodeDrive::Voltage(3.0),
            },
        ],
    };
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(300.0);
    faces[1] = ThermalFace::Temperature(300.0);
    let thermal = ThermalSetup::new(faces);
    // The domain is all solid: these unit fluid coefficients are unused by
    // the physical solve and do not purport to represent atmospheric air.
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1.0,
        kinematic_viscosity_m2_s: 1.0,
    };
    let flow = FlowField::quiescent(&domain);
    let requirement = ScalarPropertyRequirement::try_new(
        THERMAL_CONDUCTIVITY_PROPERTY,
        THERMAL_CONDUCTIVITY_DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )
    .unwrap();
    let mut electrical_states = Vec::new();
    let mut thermal_states = Vec::new();
    let result = solve_material_electrothermal(
        &domain,
        &electric,
        &thermal,
        &[300.0; 3],
        MaterialCouplingConfig {
            max_passes: 100,
            temperature_tolerance_k: 1e-9,
            power_tolerance_w: 1e-9,
            relative_tolerance: 1e-10,
        },
        |_, _, temperature| {
            let resolved = resolve_conductor_state(
                card,
                &point(temperature, 101_325.0),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .map_err(|e| source_error("material.electrical_resistivity", e))?;
            let rho = resolved.resistivity_ohm_m();
            electrical_states.push(resolved.material().clone());
            Ok(rho)
        },
        |heated, temperatures| {
            let mut conductivity = None;
            for &temperature in temperatures {
                let state = resolve_material_state_point(
                    card,
                    &point(temperature, 101_325.0),
                    core::slice::from_ref(&requirement),
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .map_err(|e| source_error("material.thermal_conductivity", e))?;
                let value = state.property(THERMAL_CONDUCTIVITY_PROPERTY).unwrap().value_si();
                if let Some(previous) = conductivity {
                    assert_eq!(previous, value, "this fixture declares constant conductivity");
                }
                conductivity = Some(value);
                thermal_states.push(state);
            }
            solve_energy(
                &domain,
                &fluid,
                &[SolidMaterial::new("resolved source coefficient", conductivity.unwrap())],
                &flow,
                heated,
                &EnergyConfig::default(),
                &gate,
            )
        },
        &gate,
    )?;
    Ok(Run { result, electrical_states, thermal_states })
}

#[test]
fn g1_g3_card_only_swaps_propagate_through_current_joule_heat_and_temperature() {
    let mut peaks = Vec::new();
    for (name, rho_scale, k) in [("fixture A", 1.0, 1.0), ("fixture B", 2.0, 2.0)] {
        let source = card(name, rho_scale, k, 320.0);
        let run = run(&source).unwrap();
        let temperature = run.result.thermal.temperature[1];
        // The two 1.5-m thermal paths are in parallel. Together with the
        // source resistivity curve, d*(1+.05*d)=3.375/(rho_scale*k).
        let (mut lo, mut hi) = (0.0, 4.0);
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if mid * (1.0 + 0.05 * mid) < 3.375 / (rho_scale * k) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        assert!((temperature - 300.0 - 0.5 * (lo + hi)).abs() < 1e-8);
        let expected_power = 9.0 / (2.0 * rho_scale * (1.0 + 0.05 * (temperature - 300.0)));
        assert!((run.result.electric.total_joule_w - expected_power).abs() < 1e-9);
        assert!((run.result.electric.delivered_w - expected_power).abs() < 1e-9);
        assert!(run.result.thermal.report.balance.relative_residual < 1e-9);
        assert!(run.result.report.joule_source_change_w < 2e-9);
        assert!(run.electrical_states.len() > 1);
        assert!(!run.thermal_states.is_empty());
        for state in run.electrical_states.iter().chain(&run.thermal_states) {
            assert_eq!(state.card_identity(), source.content_hash());
            assert!(state.query_point().iter().any(|(axis, value)| {
                axis == "pressure" && *value == 101_325.0
            }));
            assert_eq!(state.properties().len(), 1);
        }
        let final_state = run.electrical_states.last().unwrap();
        assert!(final_state.query_point().iter().any(|(axis, value)| {
            axis == "T" && *value == temperature
        }));
        assert_ne!(run.electrical_states[0].identity(), final_state.identity());
        peaks.push(temperature);
    }
    assert!(peaks[1] < peaks[0]);
}

#[test]
fn g5_spatial_and_uniform_conductors_share_the_exact_material_state() {
    let source = card("geometry-independent fixture", 1.0, 1.0, 320.0);
    let query = point(310.0, 101_325.0);
    let bulk = resolve_conductor_state(
        &source, &query, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    assert_eq!(bulk.resistivity_ohm_m(), 3.0);
    let wire = bulk.clone().into_uniform_conductor(2.0, 0.5).unwrap();
    assert_eq!(wire.resistance_ohm(), 12.0);
    assert_eq!(wire.material(), bulk.material());
    assert_eq!(wire, resolve_uniform_conductor(
        &source, &query, MaterialPropertySelection::SingleClaimOnly, 2.0, 0.5,
    ).unwrap());
    let other_pressure = resolve_conductor_state(
        &source, &point(310.0, 105_000.0), MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    assert_eq!(bulk.resistivity_ohm_m(), other_pressure.resistivity_ohm_m());
    assert_ne!(bulk.material().identity(), other_pressure.material().identity());
    for (length, area) in [(0.0, 1.0), (1.0, -1.0), (f64::NAN, 1.0), (1.0, f64::INFINITY)] {
        assert!(matches!(
            bulk.clone().into_uniform_conductor(length, area),
            Err(ConductorError::InvalidInput { .. })
        ));
    }
}

#[test]
fn g0_final_temperature_cannot_escape_the_source_curve_support() {
    let source = card("narrow-support fixture", 1.0, 1.0, 300.5);
    assert!(matches!(
        run(&source),
        Err(ChtError::InvalidInput { field: "material.electrical_resistivity", .. })
    ));
}

#[test]
fn g0_missing_or_outside_state_coordinates_are_not_replaced_by_defaults() {
    let source = card("complete-state fixture", 1.0, 1.0, 320.0);
    for query in [
        QueryPoint::new().with("T", 310.0).unwrap(),
        point(310.0, 90_000.0),
        point(330.0, 101_325.0),
    ] {
        assert!(matches!(
            resolve_conductor_state(&source, &query, MaterialPropertySelection::SingleClaimOnly),
            Err(ConductorError::MaterialState(_))
        ));
    }
}
