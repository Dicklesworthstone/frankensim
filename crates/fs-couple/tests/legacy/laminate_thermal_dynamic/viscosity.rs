//! Synthetic source → original midpoint mechanics → original finite heat.

use super::*;
use fs_couple::laminate_thermal::reciprocal::dynamic::LaminateViscosityControls;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    InterpolationPolicy, MaterialCard, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    UncertaintyModel,
};
use fs_material::plane_stress_viscosity::{
    MaterialPlane, PLANE_STRESS_VISCOSITY_PROPERTIES as KEYS, PlaneStressViscosity,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::{DynViscosity, Temperature};
use std::sync::Arc;

fn source_card(name: &str, viscosity: f64, slope_per_k: f64) -> Arc<MaterialCard> {
    let elastic = card_with_coefficients(name, 0.0, 0.0, 0.0, CP);
    let mut claims = elastic.claims().clone();
    for (index, key) in KEYS.into_iter().enumerate() {
        // Engineering gamma=2*epsilon_xy: the nu=0 isotropic axial oracle
        // requires V66=V11/2, so it is invariant under the facet-frame rotation.
        let value = match index {
            0 | 1 => viscosity,
            2 => 0.5 * viscosity,
            _ => 0.0,
        };
        let varying = slope_per_k != 0.0;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, DynViscosity::DIMS),
                value: if varying {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![
                            (280.0, value * (1.0 - 20.0 * slope_per_k)),
                            (360.0, value * (1.0 + 60.0 * slope_per_k)),
                        ],
                        dims: DynViscosity::DIMS,
                    }
                } else {
                    PropertyValue::Scalar {
                        value,
                        dims: DynViscosity::DIMS,
                    }
                },
                interpolation: if varying {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", temperature(), 280.0, 360.0)
                    .with("omega", 0.0, 20000.0),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic layered Kelvin heat oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    Arc::new(MaterialCard::assemble(elastic.id().clone(), claims, vec![]).unwrap())
}

fn viscosity(card: &MaterialCard) -> PlaneStressViscosity {
    let requirements: Vec<_> = KEYS
        .iter()
        .map(|name| {
            ScalarPropertyRequirement::try_new(
                *name,
                DynViscosity::DIMS,
                ScalarAdmissibility::Finite,
            )
            .unwrap()
        })
        .collect();
    let source = resolve_material_state_point(
        card,
        &at(T0).with("omega", 1000.0).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    PlaneStressViscosity::from_material(source, MaterialPlane::Axes12, KEYS).unwrap()
}

fn loss_controls() -> LaminateViscosityControls {
    LaminateViscosityControls {
        angular_frequency_band_rad_s: [0.0, 10000.0],
        maximum_assembly_terms: 10_000_000,
    }
}

#[test]
fn g1_kelvin_decay_returns_each_layers_original_physical_heat_once() {
    // Symmetric source loss makes B=0, preserving the existing independent
    // axial scalar oracle while requiring unequal 1:4:1 cell heat routing.
    let eta = [1e6, 4e6, 1e6];
    let cards: Vec<_> = eta
        .iter()
        .enumerate()
        .map(|(i, &v)| source_card(&format!("layer {i}"), v, 0.0))
        .collect();
    let body = prepare_body(&cards, &[1, 1, 1], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let mut thermal = initial(&bound, &[T0; 3]);
    let shell = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        axial_triangle(),
        limits(),
        || false,
    )
    .unwrap();
    let equilibrium = ReciprocalLayeredThermalShell::try_new(
        &shell,
        &thermal,
        &[0.0; 18],
        &[[280.0, 360.0]; 3],
        0.02,
        controls(),
        || false,
    )
    .unwrap();
    let mut velocity = [0.0; 18];
    velocity[6] = 0.5;
    velocity[12] = 0.5;
    let mut coupled = equilibrium
        .begin_transient(&velocity, &[[0.0; 2]], 1e-7, || false)
        .unwrap()
        .with_layer_viscosity(
            cards.iter().map(|card| viscosity(card)).collect(),
            loss_controls(),
            || false,
        )
        .unwrap();
    let volume = AREA * H;
    let mass = RHO * volume / 2.0;
    let k = E * volume;
    let c = volume * eta.iter().sum::<f64>() / 3.0;
    let dt = 2e-4;
    let h = 0.5 * dt;
    let initial_energy = thermal.energy_j() + field_internal(coupled.state());
    let mut workspace = LayeredThermalWorkspace::new(3).unwrap();
    let (mut q, mut v) = (0.0, 0.5);
    let mut total_heat = 0.0;
    for index in 1..=8 {
        let v1 = ((mass - h * c - h * h * k) * v - dt * k * q) / (mass + h * c + h * h * k);
        let mean = 0.5 * (v + v1);
        let q1 = q + dt * mean;
        let step = coupled
            .advance_to(
                &thermal,
                index as f64 * dt,
                &[0.0; 18],
                conduction(),
                LayeredThermalSolveControls::default(),
                controls(),
                &mut workspace,
                || false,
            )
            .unwrap();
        close(coupled.state().nodal_velocity_6n()[6], v1, 2e-9);
        close(coupled.state().response().displacement_6n[6], q1, 2e-11);
        let report = &step.report.thermoelastic;
        for (cell, &eta) in report.cells.iter().zip(&eta) {
            close(
                cell.dissipated_energy_j,
                dt * volume / 3.0 * eta * mean * mean,
                2e-9,
            );
            assert!(cell.dissipated_energy_j >= 0.0);
            close(cell.reversible_heat_j, 0.0, 1e-8);
            close(cell.applied_heat_j, cell.dissipated_energy_j, 1e-8);
        }
        close(report.dissipated_energy_j, dt * c * mean * mean, 5e-9);
        close(
            step.report.native.dissipated_energy_j,
            report.dissipated_energy_j,
            5e-11,
        );
        close(report.external_heat_j, 0.0, 1e-8);
        close(report.combined_energy_residual_j, 0.0, 1e-7);
        close(
            step.thermal.energy_j() + field_internal(coupled.state()),
            initial_energy,
            3e-7,
        );
        total_heat += report.dissipated_energy_j;
        thermal = step.thermal;
        (q, v) = (q1, v1);
    }
    assert!(total_heat > 0.1);
    assert!(coupled.state().temperatures_k().iter().all(|t| *t > T0));
}

#[test]
fn g1_temperature_dependent_viscosity_uses_pinned_midpoint_source_and_actual_heat() {
    let eta = 2e6;
    let slope = 0.01;
    let cards = [source_card(
        "temperature-dependent Kelvin oracle",
        eta,
        slope,
    )];
    let body = prepare_body(&cards, &[1], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let thermal = initial(&bound, &[T0]);
    let shell = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        axial_triangle(),
        limits(),
        || false,
    )
    .unwrap();
    let equilibrium = ReciprocalLayeredThermalShell::try_new(
        &shell,
        &thermal,
        &[0.0; 18],
        &[[280.0, 360.0]],
        0.02,
        controls(),
        || false,
    )
    .unwrap();
    let original = viscosity(&cards[0]);
    let selected: Vec<_> = original
        .source()
        .properties()
        .iter()
        .map(|property| property.answer().receipt.selected)
        .collect();
    let mut velocity = [0.0; 18];
    velocity[6] = 5.0;
    velocity[12] = 5.0;
    let mut coupled = equilibrium
        .begin_transient(&velocity, &[[0.0; 2]], 1e-7, || false)
        .unwrap()
        .with_layer_viscosity(vec![original], loss_controls(), || false)
        .unwrap();
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    let dt = 1e-3;
    let original_identity = coupled.layer_viscosity_sources().unwrap()[0]
        .source()
        .identity();
    let original_coordinates = coupled.state().motion().coordinates().to_vec();
    let mut one_trial = controls();
    one_trial.maximum_iterations = 1;
    assert!(matches!(
        coupled.advance_to(
            &thermal,
            dt,
            &[0.0; 18],
            conduction(),
            LayeredThermalSolveControls::default(),
            one_trial,
            &mut workspace,
            || false
        ),
        Err(LaminateThermoelasticError::NotConverged { .. })
    ));
    assert_eq!(coupled.time_s(), 0.0);
    assert_eq!(coupled.state().temperatures_k(), &[T0]);
    assert_eq!(coupled.state().motion().coordinates(), original_coordinates);
    assert_eq!(
        coupled.layer_viscosity_sources().unwrap()[0]
            .source()
            .identity(),
        original_identity
    );
    let step = coupled
        .advance_to(
            &thermal,
            dt,
            &[0.0; 18],
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut workspace,
            || false,
        )
        .unwrap();
    let mean_t = 0.5 * (T0 + coupled.state().temperatures_k()[0]);
    let prepared = &coupled.layer_viscosity_sources().unwrap()[0];
    close(
        prepared.source().state_coordinate("T").unwrap(),
        mean_t,
        1e-13,
    );
    assert_eq!(prepared.source().card_identity(), cards[0].content_hash());
    assert_eq!(
        prepared
            .source()
            .properties()
            .iter()
            .map(|property| property.answer().receipt.selected)
            .collect::<Vec<_>>(),
        selected
    );
    let midpoint_eta = eta * (1.0 + slope * (mean_t - T0));
    close(prepared.matrix_pa_s()[0], midpoint_eta, 1e-8);
    assert!(midpoint_eta > eta);
    let volume = AREA * H;
    let mass = RHO * volume / 2.0;
    let k = E * volume;
    let c = volume * midpoint_eta;
    let h = 0.5 * dt;
    let v1 = (mass - h * c - h * h * k) * 5.0 / (mass + h * c + h * h * k);
    let mean_v = 0.5 * (5.0 + v1);
    let heat = dt * c * mean_v * mean_v;
    close(coupled.state().nodal_velocity_6n()[6], v1, 1e-8);
    close(step.report.thermoelastic.dissipated_energy_j, heat, 1e-7);
    close(step.thermal.energy_j() - thermal.energy_j(), heat, 2e-7);
    close(
        RHO * volume * CP * (coupled.state().temperatures_k()[0] - T0),
        heat,
        2e-7,
    );
}

#[test]
fn g3_zero_loss_is_the_original_trajectory_and_failed_trials_leave_no_heat() {
    let cards = [source_card("zero Kelvin oracle", 0.0, 0.0)];
    let body = prepare_body(&cards, &[1], 1.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, ids()).unwrap();
    let thermal = initial(&bound, &[T0]);
    let shell = LayeredThermalShell::try_new(
        &bound,
        mechanical(&cards),
        ReferenceSurface::Midsurface,
        axial_triangle(),
        limits(),
        || false,
    )
    .unwrap();
    let mut velocity = [0.0; 18];
    velocity[6] = 0.05;
    velocity[12] = 0.05;
    let launch = || {
        ReciprocalLayeredThermalShell::try_new(
            &shell,
            &thermal,
            &[0.0; 18],
            &[[280.0, 360.0]],
            0.02,
            controls(),
            || false,
        )
        .unwrap()
        .begin_transient(&velocity, &[[0.0; 2]], 1e-7, || false)
        .unwrap()
    };
    let alien = source_card("different source", 0.0, 0.0);
    assert!(
        launch()
            .with_layer_viscosity(vec![viscosity(&alien)], loss_controls(), || false)
            .is_err()
    );
    let wrong_plane = PlaneStressViscosity::from_material(
        viscosity(&cards[0]).source().clone(),
        MaterialPlane::Axes13,
        KEYS,
    )
    .unwrap();
    assert!(
        launch()
            .with_layer_viscosity(vec![wrong_plane], loss_controls(), || false)
            .is_err()
    );
    let mut unsupported = loss_controls();
    unsupported.angular_frequency_band_rad_s[1] = 30000.0;
    assert!(
        launch()
            .with_layer_viscosity(vec![viscosity(&cards[0])], unsupported, || false)
            .is_err()
    );
    let mut bare = launch();
    let mut coupled = launch()
        .with_layer_viscosity(vec![viscosity(&cards[0])], loss_controls(), || false)
        .unwrap();
    let before = coupled.state().motion().coordinates().to_vec();
    let before_v = coupled.state().motion().velocity().to_vec();
    let receipt = coupled.layer_viscosity_sources().unwrap()[0]
        .source()
        .identity();
    let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
    let mut polls = 0;
    assert!(
        coupled
            .advance_to(
                &thermal,
                2e-4,
                &[0.0; 18],
                conduction(),
                LayeredThermalSolveControls::default(),
                controls(),
                &mut workspace,
                || {
                    polls += 1;
                    polls > 24
                }
            )
            .is_err()
    );
    assert_eq!(coupled.state().motion().coordinates(), before);
    assert_eq!(coupled.state().motion().velocity(), before_v);
    assert_eq!(coupled.state().temperatures_k(), &[T0]);
    assert_eq!(coupled.time_s(), 0.0);
    assert_eq!(
        coupled.layer_viscosity_sources().unwrap()[0]
            .source()
            .identity(),
        receipt
    );
    let zero = coupled
        .advance_to(
            &thermal,
            2e-4,
            &[0.0; 18],
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut workspace,
            || false,
        )
        .unwrap();
    let baseline = bare
        .advance_to(
            &thermal,
            2e-4,
            &[0.0; 18],
            conduction(),
            LayeredThermalSolveControls::default(),
            controls(),
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(
        coupled.state().motion().coordinates(),
        bare.state().motion().coordinates()
    );
    assert_eq!(
        coupled.state().motion().velocity(),
        bare.state().motion().velocity()
    );
    assert_eq!(
        zero.thermal.energy_j().to_bits(),
        baseline.thermal.energy_j().to_bits()
    );
    assert_eq!(zero.report.thermoelastic.dissipated_energy_j, 0.0);
}
