//! Exact source paths and reversible energy around the original J2 owner.
use super::*;
use fs_material::SmallStrainLaw;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::thermal::ThermoplasticSource;
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA,
};

fn curve_fixture(cp: f64, knots: Vec<(f64, f64)>) -> Fixture {
    let mut data = properties(cp, false);
    let mut expansion = scalar(ALPHA, INVERSE_TEMPERATURE_DIMS, 0.0);
    curve(&mut expansion, "T", knots);
    data.push(expansion);
    Fixture::from_properties(
        data,
        cp,
        QueryPoint::new().with("T", 300.0).unwrap(),
        "T",
        SolidLiquidPhase::Solid,
    )
}
fn source(fixture: &Fixture) -> ThermoplasticSource<'_> {
    let mut query = fixture.query();
    query.thermal_expansion = Some(MaterialPropertySelection::SingleClaimOnly);
    ThermoplasticSource::try_new(query, 0.0).unwrap()
}
fn energy_tolerances() -> EnergyTolerances {
    EnergyTolerances {
        energy_j_m3: 1e-5,
        yield_pa: 1e-3,
    }
}

#[test]
fn source_knots_negative_expansion_and_cooling_share_the_native_work_identity() {
    let cp = 500.0;
    let fixture = curve_fixture(
        cp,
        vec![
            (280.0, -1e-5),
            (300.0, 1e-5),
            (306.0, 3e-5),
            (312.0, 1e-5),
            (320.0, -1e-5),
        ],
    );
    let source = source(&fixture);
    let virgin = source.reference().initial_state();
    // Three independent trapezoids: 6*(1+3)e-5/2 +
    // 6*(3+1)e-5/2 + 4*(1+0)e-5/2. The endpoint alpha is ZERO.
    let free = 2.6e-4;
    for component in &source.mechanical_strain(316.0, &[0.0; 6]).unwrap()[..3] {
        close(*component, -free, 1e-18);
    }
    let heated = source
        .advance_with_energy(
            0.0,
            cp * 16.0,
            &[0.0; 6],
            &virgin,
            &[0.0; 6],
            energy_tolerances(),
        )
        .unwrap();
    let nine_bulk = 6e9;
    let trace_sigma = -nine_bulk * free;
    for component in &heated.mechanical.stress_pa[..3] {
        close(*component, trace_sigma / 3.0, 1e-8);
    }
    // s1=alpha(316)*tr(sigma1)=0, so all reversible heat is
    // tr(sigma1)*delta(f), even though the integrated free strain is nonzero.
    close(heated.thermoelastic_storage_change_j_m3, 0.0, 1e-8);
    close(heated.reversible_heat_j_m3, trace_sigma * free, 1e-8);
    close(
        heated.mechanical.energy.physical_dissipation_j_m3,
        0.0,
        1e-12,
    );
    close(
        heated.mechanical.energy.final_state.elastic_j_m3,
        0.5 * nine_bulk * free * free,
        1e-8,
    );
    close(
        heated.mechanical.energy.numerical_dissipation_j_m3,
        0.5 * nine_bulk * free * free,
        1e-8,
    );

    let cooled = source
        .advance_with_energy(
            cp * 16.0,
            0.0,
            &[0.0; 6],
            &heated.mechanical.state,
            &[0.0; 6],
            energy_tolerances(),
        )
        .unwrap();
    assert_eq!(cooled.mechanical.state, virgin);
    for stress in cooled.mechanical.stress_pa {
        close(stress, 0.0, 1e-8);
    }
    close(cooled.mechanical.energy.final_state.total_j_m3(), 0.0, 1e-8);
    // A finite backward-Euler cycle retains its native numerical loss;
    // it cannot be renamed plastic heat to manufacture reversibility.
    close(
        heated.reversible_heat_j_m3
            + cooled.reversible_heat_j_m3
            + heated.mechanical.energy.numerical_dissipation_j_m3
            + cooled.mechanical.energy.numerical_dissipation_j_m3,
        0.0,
        1e-8,
    );

    // Actual free expansion and its reverse have zero stress/storage/heat.
    let free_strain = [free, free, free, 0.0, 0.0, 0.0];
    let free_heating = source
        .advance_with_energy(
            0.0,
            cp * 16.0,
            &[0.0; 6],
            &virgin,
            &free_strain,
            energy_tolerances(),
        )
        .unwrap();
    let free_cooling = source
        .advance_with_energy(
            cp * 16.0,
            0.0,
            &free_strain,
            &free_heating.mechanical.state,
            &[0.0; 6],
            energy_tolerances(),
        )
        .unwrap();
    for increment in [&free_heating, &free_cooling] {
        for stress in increment.mechanical.stress_pa {
            close(stress, 0.0, 1e-8);
        }
        close(increment.caloric_heat_j_m3(), 0.0, 1e-8);
        close(
            increment.mechanical.energy.numerical_dissipation_j_m3,
            0.0,
            1e-8,
        );
    }
    assert_eq!(free_cooling.mechanical.state, virgin);
}

#[test]
fn a_surface_path_uses_its_crossed_cells_and_explicit_temperature_axis() {
    let axis = "specimen_kelvin";
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let mut data = properties(500.0, false);
    let mut expansion = scalar(ALPHA, INVERSE_TEMPERATURE_DIMS, 0.0);
    pressure_surface(&mut expansion, [1e-5, 4e-5, 2e-5, 6e-5]);
    data.push(expansion);
    for property in &mut data {
        match &mut property.value {
            PropertyValue::Curve { abscissa, .. } => {
                *abscissa = axis.into();
                property.validity = ValidityDomain::unconstrained().with(axis, 280.0, 320.0);
            }
            PropertyValue::Surface2D { axes, .. } => {
                axes[0] = axis.into();
                property.validity = ValidityDomain::unconstrained()
                    .with(axis, 280.0, 320.0)
                    .with("pressure", 1e5, 3e5);
            }
            _ => {}
        }
    }
    let point_context = QueryPoint::new()
        .with_quantity(axis, absolute, 300.0)
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            2e5,
        )
        .unwrap();
    let fixture =
        Fixture::from_properties(data, 500.0, point_context, axis, SolidLiquidPhase::Solid);
    let source = source(&fixture);
    // At p=2e5: alpha is 1.5e-5 at 280 K, 3e-5 at the triangle
    // boundary 300 K, then 5e-5 at 320 K. Both one-sided slopes matter.
    let hot = source.mechanical_strain(310.0, &[0.0; 6]).unwrap();
    let cold = source.mechanical_strain(290.0, &[0.0; 6]).unwrap();
    close(hot[0], -3.5e-4, 1e-18);
    close(cold[0], 2.625e-4, 1e-18);
    let (_, resolved) = source.resolve_enthalpy(5000.0).unwrap();
    let reference = source.reference().material().property(ALPHA).unwrap();
    let current = resolved.material().property(ALPHA).unwrap();
    close(current.value_si(), 4e-5, 1e-18);
    assert_eq!(
        current.answer().receipt.selected,
        reference.answer().receipt.selected
    );
    assert_eq!(
        resolved.material().axis_quantities().get(axis),
        Some(&absolute)
    );
    assert!(
        !resolved
            .material()
            .query_point()
            .iter()
            .any(|(name, _)| name == "T")
    );
    let energy = source
        .advance_with_energy(
            0.0,
            5000.0,
            &[0.0; 6],
            &source.reference().initial_state(),
            &[0.0; 6],
            energy_tolerances(),
        )
        .unwrap();
    let trace_sigma = -6e9 * 3.5e-4;
    let storage = 310.0 * 4e-5 * trace_sigma;
    close(energy.thermoelastic_storage_change_j_m3, storage, 1e-7);
    close(
        energy.reversible_heat_j_m3,
        trace_sigma * 3.5e-4 - storage,
        1e-7,
    );
}

#[test]
fn an_interior_alpha_slope_cannot_bypass_capacity_admission_at_a_zero_reference() {
    let fixture = curve_fixture(
        500.0,
        vec![
            (280.0, 0.0),
            (300.0, 0.0),
            (301.0, 1e-4),
            (302.0, 0.0),
            (320.0, 0.0),
        ],
    );
    let original = fixture.caloric.curve().clone();
    let mut query = fixture.query();
    query.thermal_expansion = Some(MaterialPropertySelection::SingleClaimOnly);
    assert!(matches!(
        CaloricPoint::try_new(query, 0.0, VOLUME, SURFACE),
        Err(FeedbackError::Invalid(message)) if message.contains("effective heat capacity")
    ));
    assert_eq!(fixture.caloric.curve(), &original);

    // This is a sufficient bound on the SAME norm admitted by the energy owner.
    // A hydrostatic elastic strain saturates |tr(e)| <= sqrt(3)||e||.
    let limit = 0.04;
    let elastic = [
        limit / 3.0_f64.sqrt(),
        limit / 3.0_f64.sqrt(),
        limit / 3.0_f64.sqrt(),
        0.0,
        0.0,
        0.0,
    ];
    close(
        fs_material::tensor::contract(&elastic, &elastic),
        limit * limit,
        1e-18,
    );
    let bulk = 1.2e9 / (3.0 * (1.0 - 0.4));
    let trace_sigma = 3.0 * bulk * fs_material::tensor::trace(&elastic);
    close(trace_sigma, 3.0 * 3.0_f64.sqrt() * bulk * limit, 1e-7);
    let effective_capacity =
        1000.0 * 500.0 + 301.5 * (-1e-4) * trace_sigma - 9.0 * bulk * 301.5 * (5e-5_f64).powi(2);
    assert!(effective_capacity < 0.0);
}
