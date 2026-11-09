//! Explicit liquid wall-Pr correction on the existing cell and link heat owners.
use super::*;

#[path = "liquid_wall/conjugate.rs"]
mod conjugate;
use fs_material::liquid::{LiquidPropertyRequirements, ResolvedLiquidState};

fn corrected_law() -> SourceCellConvectionLaw {
    SourceCellConvectionLaw::ForcedTubeSectionLiquidWall {
        speed_m_s: 1.0,
        start: Length::new(1.0),
        end: Length::new(2.0),
    }
}

fn corrected_link() -> SourceConvectiveLink {
    SourceConvectiveLink {
        property_temperature: ConvectionPropertyTemperature::Reservoir,
        law: corrected_law(),
        ..link()
    }
}

// Independent Gnielinski Eq. (2), corrected friction relation and Eq. (3).
// The moderate interval permits ordinary differences of cumulative transfer.
// Fixture Cp is constant; k and mu are the original temperature curves.
fn coefficient(bulk: f64, wall: f64, corrected: bool) -> f64 {
    let mu = 0.001 * (1.0 - 0.01 * (bulk - 300.0));
    let k = 0.5 + 0.002 * (bulk - 300.0);
    let re = 1000.0 * 0.1 / mu;
    let pr = mu * FLUID_CP / k;
    let wall_pr = 0.001 * (1.0 - 0.01 * (wall - 300.0)) * FLUID_CP
        / (0.5 + 0.002 * (wall - 300.0));
    let f = (1.8 * re.log10() - 1.5).powi(-2);
    let developed = (f / 8.0) * (re - 1000.0) * pr
        / (1.0 + 12.7 * (f / 8.0).sqrt() * (pr.powf(2.0 / 3.0) - 1.0));
    let section = ((20.0 + 20.0_f64.cbrt()) - (10.0 + 10.0_f64.cbrt())) / 10.0;
    let correction = if corrected { (pr / wall_pr).powf(0.11) } else { 1.0 };
    developed * section * correction * k / 0.1
}

#[test]
fn liquid_wall_pr_changes_heating_and_cooling_using_original_bulk_and_wall_sources() {
    let chart = curve();
    let thermal = initial(&chart, [310.0, 290.0], 0.0, 0.0);
    let (card, _) = fluid_source(1.0);
    let liquid = liquid_law(&card);
    let bindings = IDS.map(|cell| SourceConvectiveCell {
        cell,
        reference_temperature: Temperature::new(300.0),
        characteristic_length: Length::new(0.1),
        property_temperature: ConvectionPropertyTemperature::Reservoir,
        law: corrected_law(),
    });
    let mut queries = Vec::new();
    let corrected = evaluate_source_cell_convection(
        &thermal, &bindings, 2,
        |binding, temperature| {
            queries.push((binding.cell, temperature.value()));
            liquid.at_temperature(temperature)
        }, || false,
    ).unwrap();
    assert_eq!(queries.len(), 4);
    for sample in corrected.convection().samples() {
        let wall = sample.wall_temperature_k;
        let expected = coefficient(300.0, wall, true);
        near(sample.heat_into_cell_w, AREA * expected * (300.0 - wall), 2e-10);
        let evaluation = sample.evaluation.as_ref().unwrap();
        near(evaluation.coefficient.qoi, expected, 2e-9);
        near(evaluation.effective_nusselt(), expected * 0.1 / 0.5, 4e-10);
        let section = evaluation.tube_section.as_ref().unwrap();
        assert_eq!((section.start.value(), section.end.value()), (1.0, 2.0));
        assert!(section.upstream.is_some());
        assert!(evaluation.coefficient.model.discrepancy_rel.is_infinite());
        let retained = corrected.materials().iter().find(|m| m.cell == sample.cell).unwrap();
        let surface = retained.surface.as_ref().unwrap();
        assert_eq!(retained.property_temperature_k, 300.0);
        assert_eq!(surface.thermal.temperature_k(), wall);
        assert_eq!(retained.material.material().card_identity(), card.content_hash());
        assert_eq!(surface.material.material().card_identity(), card.content_hash());
        assert!(queries.contains(&(sample.cell, 300.0)));
        assert!(queries.contains(&(sample.cell, wall)));
        let borrowed = &retained.material;
        assert!(core::ptr::eq(
            ThermalFluidTransport::resolved_liquid_source(&borrowed).unwrap(),
            &retained.material,
        ));
    }
    let uncorrected = bindings.map(|binding| SourceConvectiveCell {
        law: SourceCellConvectionLaw::ForcedTubeSection {
            speed_m_s: 1.0, start: Length::new(1.0), end: Length::new(2.0),
        },
        ..binding
    });
    let mut calls = 0;
    let held = evaluate_source_cell_convection(&thermal, &uncorrected, 2, |_, t| {
        calls += 1;
        liquid.at_temperature(t)
    }, || false).unwrap();
    assert_eq!(calls, 2, "unselected K=1 keeps its single bulk query");
    assert!(held.materials().iter().all(|sample| sample.surface.is_none()));
    for sample in held.convection().samples() {
        near(sample.evaluation.as_ref().unwrap().coefficient.qoi,
            coefficient(300.0, sample.wall_temperature_k, false), 2e-9);
    }
    assert!(corrected.convection().sources_w()[&IDS[0]].abs()
        > held.convection().sources_w()[&IDS[0]].abs());
    assert!(corrected.convection().sources_w()[&IDS[1]].abs()
        < held.convection().sources_w()[&IDS[1]].abs());
}

fn closed_endpoint(dt: f64, corrected: bool) -> [f64; 2] {
    let (mut lo, mut hi) = (294.0, 310.0);
    for _ in 0..80 {
        let wall = 0.5 * lo + 0.5 * hi;
        let bulk = 290.0 + CP / FLUID_CP * (310.0 - wall);
        let heat = dt * LINK_AREA * coefficient(bulk, wall, corrected) * (wall - bulk);
        if MASS * CP * (wall - 310.0) + heat > 0.0 { hi = wall; }
        else { lo = wall; }
    }
    let wall = 0.5 * lo + 0.5 * hi;
    [wall, 290.0 + CP / FLUID_CP * (310.0 - wall)]
}

#[test]
fn liquid_wall_pr_drives_both_original_finite_stores_without_external_heat() {
    let chart = curve();
    let (card, caloric) = fluid_source(1.0);
    let liquid = liquid_law(&card);
    let dt = 0.2;
    let expected = closed_endpoint(dt, true);
    let held = closed_endpoint(dt, false);
    assert!(expected[0] < held[0] - 0.01);
    assert!(expected[1] > held[1] + 0.002);
    for reversed in [false, true] {
        let thermal = pair(&chart, &caloric, reversed, 0.0);
        let saved = thermal.cells().clone();
        let step = BackwardEulerStep::new(&thermal, dt, || Ok(())).unwrap();
        let increments = step.cell_ids().iter().map(|&id| {
            if id == IDS[0] { MASS * CP * (expected[0] - 310.0) }
            else { MASS * FLUID_CP * (expected[1] - 290.0) }
        }).collect::<Vec<_>>();
        let mut queries = Vec::new();
        let trial = step.evaluate(&increments, |state| {
            let resolved = evaluate_source_link_convection(
                state, &[corrected_link()], 1,
                |_, _, temperature| {
                    queries.push(temperature.value());
                    liquid.at_temperature(temperature)
                }, || false,
            )?;
            Ok::<_, SourceCellConvectionError<MaterialStatePointError>>(StepCoefficients {
                sources_w: BTreeMap::new(),
                conductances_w_k: resolved.conductances_w_k().clone(),
                payload: resolved,
            })
        }, || Ok(())).unwrap();
        let (next, resolved, audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
        assert_eq!(queries.len(), 2);
        near(queries[0], expected[1], 1e-12);
        near(queries[1], expected[0], 1e-12);
        let sample = &resolved.samples()[0];
        near(sample.property_temperature_k, expected[1], 1e-12);
        near(sample.surface.as_ref().unwrap().thermal.temperature_k(), expected[0], 1e-12);
        near(sample.conductance_w_k, LINK_AREA * coefficient(expected[1], expected[0], true), 1e-10);
        near(MASS * CP * (expected[0] - 310.0) + dt * sample.heat_into_fluid_w, 0.0, 1e-8);
        near(MASS * FLUID_CP * (expected[1] - 290.0) - dt * sample.heat_into_fluid_w, 0.0, 1e-8);
        for (cell, temperature) in IDS.into_iter().zip(expected) {
            near(next.temperature_k(cell).unwrap(), temperature, 1e-9);
        }
        near(next.energy_j(), thermal.energy_j(), 1e-8);
        near(audit.external_heat_j, 0.0, 1e-12);
        near(audit.residual_j, 0.0, 1e-8);
        assert_eq!(thermal.cells(), &saved);
    }
}

struct IncoherentLiquidView {
    original: ResolvedLiquidState,
    thermal: ThermalFluidTransportState,
}
impl FluidTransport for IncoherentLiquidView {
    fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
        Ok(self.thermal.transport())
    }
}
impl ThermalFluidTransport for IncoherentLiquidView {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, fs_material::MaterialError> {
        Ok(self.thermal)
    }
    fn resolved_liquid_source(&self) -> Option<&ResolvedLiquidState> {
        Some(&self.original)
    }
}

#[test]
fn liquid_wall_pr_refuses_unproved_liquids_and_incoherent_source_witnesses() {
    let chart = curve();
    let (card, caloric) = fluid_source(1.0);
    let liquid = liquid_law(&card);
    let thermal = pair(&chart, &caloric, false, 0.0);
    let saved = thermal.cells().clone();
    let raw = evaluate_source_link_convection(&thermal, &[corrected_link()], 1,
        |_, _, t| liquid.at_temperature(t).map(|state| state.thermal_transport_state().unwrap()),
        || false,
    );
    assert!(matches!(raw, Err(SourceCellConvectionError::Convection(
        ConvectionSourceError::Invalid("liquid wall-Pr correction requires an original resolved liquid source")
    ))));
    let gas = GasSpec::dry_air_ussa1976();
    let gas_result = evaluate_source_link_convection(&thermal, &[corrected_link()], 1,
        |_, _, t| GasState::try_new(&gas, t.value(), 100_000.0), || false,
    );
    assert!(matches!(gas_result, Err(SourceCellConvectionError::Convection(
        ConvectionSourceError::Invalid("liquid wall-Pr correction requires an original resolved liquid source")
    ))));
    let incoherent = evaluate_source_link_convection(&thermal, &[corrected_link()], 1,
        |_, _, t| liquid.at_temperature(t).map(|original| {
            let state = original.transport_state().unwrap();
            let changed = FluidTransportState::try_new(
                state.density_kg_m3(), state.specific_heat_j_kg_k(),
                1.1 * state.thermal_conductivity_w_m_k(), state.dynamic_viscosity_pa_s(),
            ).unwrap();
            IncoherentLiquidView {
                original,
                thermal: ThermalFluidTransportState::try_new(changed, t.value(), None).unwrap(),
            }
        }), || false,
    );
    assert!(matches!(incoherent, Err(SourceCellConvectionError::Convection(
        ConvectionSourceError::Invalid("liquid source witness does not match the retained thermal transport state")
    ))));
    use fs_airflow::graph::thermal::medium::convection::implicit::{
        ExteriorConvectionSettings, ImplicitLiquidConvection,
    };
    let settings = ExteriorConvectionSettings {
        fluid_temperature_bounds_k: [280.0, 320.0],
        surface_temperature_bounds_k: [280.0, 320.0],
        characteristic_length: Length::new(0.1),
        property_temperature: ConvectionPropertyTemperature::Reservoir,
        law: corrected_law(),
        emissivity: 0.0,
    };
    assert!(ImplicitLiquidConvection::try_new(&chart, [-20_000.0, 20_000.0], &liquid, settings).is_err());
    let held = ExteriorConvectionSettings {
        law: SourceCellConvectionLaw::ForcedTubeSection {
            speed_m_s: 1.0, start: Length::new(1.0), end: Length::new(2.0),
        },
        ..settings
    };
    assert!(ImplicitLiquidConvection::try_new(&chart, [-20_000.0, 20_000.0], &liquid, held)
        .unwrap().with_gnielinski_wall_correction().is_ok());
    assert_eq!(thermal.cells(), &saved);
}

#[test]
fn liquid_wall_pr_keeps_selected_claims_and_fixed_non_temperature_context() {
    let chart = curve();
    let (original, caloric) = fluid_source(1.0);
    let mut claims = ClaimSet::new();
    let mut first = Vec::new();
    let mut second = Vec::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        claim.validity = claim.validity.with("p", 100_000.0, 200_000.0);
        let selected = claims.insert_claim(claim.clone()).unwrap();
        first.push((claim.key.name().to_owned(), selected));
        if claim.key.name() == SPECIFIC_HEAT_CAPACITY_PROPERTY {
            claim.value = PropertyValue::Scalar { value: 5000.0, dims: SPECIFIC_HEAT_CAPACITY_DIMS };
            second.push((claim.key.name().to_owned(), claims.insert_claim(claim).unwrap()));
        } else {
            second.push((claim.key.name().to_owned(), selected));
        }
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let make = |pressure, selection| LiquidTemperatureLaw::try_new(
        &card, &context(300.0).with("p", pressure).unwrap(), "T",
        LiquidViscosityBasis::Dynamic, MaterialPropertySelection::PinnedByProperty(selection),
    ).unwrap();
    let bulk = make(100_000.0, first.clone());
    let other_claim = make(100_000.0, second);
    let other_pressure = make(200_000.0, first);
    assert_eq!(bulk.reference_state().requirements(),
        &LiquidPropertyRequirements::canonical(LiquidViscosityBasis::Dynamic).unwrap());
    let thermal = pair(&chart, &caloric, false, 0.0);
    let saved = thermal.cells().clone();
    for changed in [&other_claim, &other_pressure] {
        let result = evaluate_source_link_convection(&thermal, &[corrected_link()], 1,
            |_, _, t| if t.value() < 300.0 {
                bulk.at_temperature(t)
            } else {
                changed.at_temperature(t)
            },
            || false,
        );
        assert!(matches!(result, Err(SourceCellConvectionError::Convection(
            ConvectionSourceError::Invalid("liquid wall-Pr correction must retain one card, selected source bundle and non-temperature context")
        ))));
    }
    assert_eq!(thermal.cells(), &saved);
}
