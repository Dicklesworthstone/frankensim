//! The original density source and reference mass survive complete melt release.
use super::*;
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;
use fs_material::fluid_transport::FluidTransport;
use fs_material::phase::transport::{
    PhaseTransportSourceQuery, resolve_phase_transport_state_with_density_departure,
};
use fs_material::thermal::source::ThermalMaterialPoint;

fn thermal_density() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, 1.0),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            2.0,
        ),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, 1.0),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            0.1,
        ),
    ] {
        let mut claim = scalar(name, dims, value);
        let liquid_only =
            name == SPECIFIC_HEAT_CAPACITY_PROPERTY || name == DYNAMIC_VISCOSITY_PROPERTY;
        claim.validity = claim.validity.with(
            "liquid",
            if liquid_only { f64::from_bits(1) } else { 0.0 },
            1.0,
        );
        if name == DENSITY_PROPERTY {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(450.0, 1.01), (500.0, 1.0), (600.0, 0.98)],
                dims,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    card(claims)
}

#[test]
fn complete_melting_keeps_source_buoyancy_and_further_heating_changes_native_motion() {
    let card = thermal_density();
    let curve = chart(&card);
    let mushy = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = domain();
    let ids = ids(&domain);
    let config = config();
    let gate = CancelGate::new_clock_free();
    let context = point_liquid(500.0);
    let liquid = LiquidTemperatureLaw::try_new(
        &card,
        &context,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let surface = resolve_liquid_surface_state(
        &card,
        &context,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let gravity = [0.2, 0.0, 0.0];
    let bounds = BoussinesqConfig {
        gravity_m_s2: gravity,
        max_relative_density_departure: 0.03,
        max_relative_specific_heat_departure: 0.0,
    };
    let original = initial(&domain, &curve, -100.0, 2000.0, 0.15);
    let history = MaterialFlowFeedback::try_new_phase_boussinesq(
        &original,
        &domain,
        FlowEnthalpyInput {
            cell_ids: &ids,
            material_id: MAT,
            caloric: &curve,
            inlet_specific_h_j_kg: [None; 6],
            max_cells: 64,
            max_connections: 512,
        },
        &config,
        (),
        &gate,
        100,
        &mushy,
        liquid.reference_state().transport_state().unwrap(),
        bounds,
        |_, cell: &DepositedCell, _| {
            resolve_phase_transport_state_with_density_departure(
                &mushy,
                cell.specific_h_j_kg,
                PhaseTransportSourceQuery {
                    source: ThermalMaterialPoint {
                        material: &card,
                        context: &context,
                        temperature_axis: "T",
                        liquid_fraction_axis: Some("liquid"),
                    },
                    viscosity_basis: LiquidViscosityBasis::Dynamic,
                    extra_requirements: &[],
                    selection: MaterialPropertySelection::SingleClaimOnly,
                },
                bounds.max_relative_density_departure,
            )
        },
    )
    .unwrap()
    .with_unsteady_viscous_heating(viscous())
    .unwrap();
    let (molten, history, _) = original
        .advance_observed(&history, 0.15, step(), || Ok(()))
        .unwrap();
    close(molten.energy_j() - original.energy_j(), 8.0 * 300.0, 1e-8);
    let before = (
        molten.cells().clone(),
        history.unsteady_flow().unwrap().clone(),
    );
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, 1.0, 0.01, 0.2).unwrap();
        assert!(
            history
                .release_fully_liquid(
                    &molten,
                    scaling,
                    receiving([0.1, 0.0, 0.0]),
                    &liquid,
                    &curve,
                    StationaryConstraintRelease::Withdraw,
                    limits(),
                    FreeSurfaceEnthalpyControls::default(),
                    cx,
                )
                .is_err()
        );
        let (released, report) = history
            .release_fully_liquid(
                &molten,
                scaling,
                receiving(gravity),
                &liquid,
                &curve,
                StationaryConstraintRelease::Withdraw,
                limits(),
                FreeSurfaceEnthalpyControls::default(),
                cx,
            )
            .unwrap();
        close(report.source.mass_kg, report.target.mass_kg, 1e-12);
        close(report.source_enthalpy_j, report.target_enthalpy_j, 1e-9);
        close(released.physical_time_s(), 0.15, 1e-15);
        close(
            released.material_at([1, 1, 0]).unwrap().density_kg_m3(),
            0.99,
            1e-12,
        );
        let mut reference = released.clone();
        let mut heated = released;
        reference.step(&gate).unwrap();
        let heat: Vec<_> = ids
            .iter()
            .enumerate()
            .filter_map(|(voxel, id)| {
                id.map(|_| FreeSurfaceCellHeat {
                    coordinates: domain.coords(voxel),
                    energy_j: 20.0,
                })
            })
            .collect();
        let accepted = heated.step_with_heat(&heat, &gate).unwrap();
        close(accepted.mass_after_kg, 8.0, 1e-11);
        close(
            accepted.enthalpy_after_j - accepted.enthalpy_before_j,
            160.0,
            1e-7,
        );
        close(
            heated.material_at([1, 1, 0]).unwrap().density_kg_m3(),
            0.988,
            1e-10,
        );
        let ordinary_speed = reference
            .flow()
            .sample_at_si(1, 1, 0)
            .unwrap()
            .unwrap()
            .velocity_m_s[0];
        let heated_speed = heated
            .flow()
            .sample_at_si(1, 1, 0)
            .unwrap()
            .unwrap()
            .velocity_m_s[0];
        assert!(heated_speed > 0.0 && heated_speed < ordinary_speed);
    });
    assert_eq!(molten.cells(), &before.0);
    assert_eq!(history.unsteady_flow().unwrap(), &before.1);
}
