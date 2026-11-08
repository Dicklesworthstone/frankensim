use super::*;
#[path = "tests/density.rs"]
mod density;
use crate::thermal::transport::{SPECIFIC_HEAT, isobaric::IsobaricHeatCapacityQuery};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyValue, Provenance,
    SelectionPolicy, UncertaintyModel,
};
use fs_qty::{Pressure, Temperature};
use fs_thermal::crystallization::CrystallizationSolve;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, TimedEvent,
    crystallization::{
        CrystallizationControls, CrystallizationMaterial, CrystallizationSimulation,
    },
};

const SELECTION: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

fn point(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
        .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), p)
        .unwrap()
}

struct Fixture {
    card: MaterialCard,
    parent_key: PropertyKey,
    crystal_key: PropertyKey,
    gap_key: PropertyKey,
    reference: QueryPoint,
}

impl Fixture {
    fn new(
        parent: [f64; 2],
        crystal: [f64; 2],
        reference_k: f64,
        gap: f64,
        gap_support: [bool; 2],
    ) -> Self {
        let parent_key = PropertyKey::new("parent_cp", SPECIFIC_HEAT);
        let crystal_key = PropertyKey::new("crystal_cp", SPECIFIC_HEAT);
        let gap_key = PropertyKey::new("parent_minus_crystal_enthalpy", SPECIFIC_ENERGY);
        let temperature = QuantitySpec::dimensional(Temperature::DIMS);
        let pressure = QuantitySpec::dimensional(Pressure::DIMS);
        let provenance = || Provenance {
            source: "synthetic phase calorimetry; numerical oracle only".into(),
            license: "test-only".into(),
            artifact: None,
        };
        let mut claims = ClaimSet::new();
        for (key, cp) in [(&parent_key, parent), (&crystal_key, crystal)] {
            claims
                .insert_claim(PropertyClaim {
                    key: key.clone(),
                    value: PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(300.0, cp[0]), (500.0, cp[1])],
                        dims: SPECIFIC_HEAT,
                    },
                    validity: ValidityDomain::unconstrained()
                        .with_quantity("T", temperature, 300.0, 500.0)
                        .with_quantity("p", pressure, 100000.0, 100000.0),
                    interpolation: InterpolationPolicy::LinearInside,
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: provenance(),
                })
                .unwrap();
        }
        let mut validity = ValidityDomain::unconstrained();
        if gap_support[0] {
            validity = validity.with_quantity("T", temperature, reference_k, reference_k);
        }
        if gap_support[1] {
            validity = validity.with_quantity("p", pressure, 100000.0, 100000.0);
        }
        claims
            .insert_claim(PropertyClaim {
                key: gap_key.clone(),
                value: PropertyValue::Scalar {
                    value: gap,
                    dims: SPECIFIC_ENERGY,
                },
                validity,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "one synthetic composition".into(),
                phase: "parent and crystal".into(),
                process: "fixed-pressure calorimetry fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        Self {
            card,
            parent_key,
            crystal_key,
            gap_key,
            reference: point(reference_k, 100000.0),
        }
    }

    fn standard() -> Self {
        // Cp_parent=3T, Cp_crystal=2T; L(400)=80000 => L(T)=T²/2.
        Self::new(
            [900.0, 1500.0],
            [600.0, 1000.0],
            400.0,
            80000.0,
            [true, true],
        )
    }

    fn caloric(&self, key: &PropertyKey, datum: f64) -> IsobaricHeatCapacityLaw {
        IsobaricHeatCapacityLaw::try_new(
            &self.card,
            IsobaricHeatCapacityQuery {
                specific_heat_cp: key,
                temperature_axis: "T",
                pressure_axis: "p",
                lower: &point(300.0, 100000.0),
                upper: &point(500.0, 100000.0),
                selection: SELECTION,
            },
            datum,
        )
        .unwrap()
    }

    fn pair<'a>(
        &'a self,
        parent: &'a IsobaricHeatCapacityLaw,
        crystal: &'a IsobaricHeatCapacityLaw,
    ) -> Result<IsobaricCrystallizationPair<'a>, ThermalTransportError> {
        IsobaricCrystallizationPair::try_new(
            &self.card,
            IsobaricCrystallizationQuery {
                parent,
                crystal,
                enthalpy_gap: &self.gap_key,
                reference: &self.reference,
                selection: SELECTION,
                temperature_range_k: [300.0, 500.0],
            },
        )
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g1_original_cp_drives_temperature_dependent_latent_energy_and_phase_inverses() {
    let fixture = Fixture::standard();
    let parent = fixture.caloric(&fixture.parent_key, 300.0);
    let crystal = fixture.caloric(&fixture.crystal_key, 450.0);
    let pair = fixture.pair(&parent, &crystal).unwrap();
    assert!(core::ptr::eq(pair.parent().source(), &parent));
    assert!(core::ptr::eq(pair.crystal().source(), &crystal));
    assert_eq!(pair.enthalpy_gap_source().evidence.value.value, 80000.0);
    for t in [300.0, 325.0, 400.0, 475.0, 500.0] {
        let parent_h = 1.5 * (t * t - 160000.0);
        let crystal_h = t * t - 240000.0;
        close(
            pair.parent().enthalpy_at_temperature(t).unwrap(),
            parent_h,
            2e-10,
        );
        close(
            pair.crystal().enthalpy_at_temperature(t).unwrap(),
            crystal_h,
            2e-10,
        );
        close(pair.latent_enthalpy_j_kg(t).unwrap(), 0.5 * t * t, 2e-10);
        close(pair.parent().temperature_k(parent_h).unwrap(), t, 2e-12);
        close(pair.crystal().temperature_k(crystal_h).unwrap(), t, 2e-12);
    }
    for branch in [pair.parent(), pair.crystal()] {
        assert_eq!(branch.storage_basis(), ThermalStorageBasis::Enthalpy);
        assert_eq!(
            CrystallizationCaloric::density_kg_m3(branch, 400.0).unwrap(),
            None
        );
        assert_eq!(
            ThermalLaw::caloric_identity(branch),
            CrystallizationCaloric::caloric_identity(branch)
        );
        let bounds = branch.specific_enthalpy_span_j_kg();
        assert!(branch.temperature_k(bounds[0].next_down()).is_err());
        assert!(branch.temperature_k(bounds[1].next_up()).is_err());
    }
    // Relative phase entropy is deliberately not inferred; these are within-branch changes.
    let a = pair.parent().enthalpy_at_temperature(325.0).unwrap();
    let b = pair.parent().enthalpy_at_temperature(475.0).unwrap();
    close(
        pair.parent()
            .specific_entropy_change_j_kg_k(a, b)
            .unwrap()
            .unwrap(),
        450.0,
        2e-10,
    );

    // The original source's arbitrary sensible zero does not change physical rebasing.
    let shifted_parent = fixture.caloric(&fixture.parent_key, 475.0);
    let shifted_crystal = fixture.caloric(&fixture.crystal_key, 325.0);
    let shifted = fixture.pair(&shifted_parent, &shifted_crystal).unwrap();
    assert_ne!(pair.identity(), shifted.identity());
    for t in [325.0, 475.0] {
        close(
            pair.parent().enthalpy_at_temperature(t).unwrap(),
            shifted.parent().enthalpy_at_temperature(t).unwrap(),
            0.0,
        );
        close(
            pair.crystal().enthalpy_at_temperature(t).unwrap(),
            shifted.crystal().enthalpy_at_temperature(t).unwrap(),
            0.0,
        );
    }
}

#[test]
fn g0_quadratic_gap_minimum_cannot_hide_between_positive_endpoint_values() {
    // Cp_parent-Cp_crystal=-1000+10*(T-300). Positive L at 300 and
    // 500 K hides a negative minimum at 400 K; neither source has a knot there.
    let bad = Fixture::new([1000.0, 3000.0], [2000.0; 2], 300.0, 10000.0, [true, true]);
    let parent = bad.caloric(&bad.parent_key, 300.0);
    let crystal = bad.caloric(&bad.crystal_key, 300.0);
    assert!(bad.pair(&parent, &crystal).is_err());
    let good = Fixture::new([1000.0, 3000.0], [2000.0; 2], 300.0, 60000.0, [true, true]);
    let parent = good.caloric(&good.parent_key, 300.0);
    let crystal = good.caloric(&good.crystal_key, 300.0);
    close(
        good.pair(&parent, &crystal)
            .unwrap()
            .latent_enthalpy_j_kg(400.0)
            .unwrap(),
        10000.0,
        1e-10,
    );
}

#[test]
fn g0_phase_sources_require_one_context_and_explicit_gap_support() {
    for support in [[false, true], [true, false]] {
        let fixture = Fixture::new([900.0, 1500.0], [600.0, 1000.0], 400.0, 80000.0, support);
        let parent = fixture.caloric(&fixture.parent_key, 300.0);
        let crystal = fixture.caloric(&fixture.crystal_key, 300.0);
        assert!(fixture.pair(&parent, &crystal).is_err());
    }
    let fixture = Fixture::standard();
    let parent = fixture.caloric(&fixture.parent_key, 300.0);
    let crystal = fixture.caloric(&fixture.crystal_key, 300.0);
    let query = IsobaricCrystallizationQuery {
        parent: &parent,
        crystal: &crystal,
        enthalpy_gap: &fixture.gap_key,
        reference: &fixture.reference,
        selection: SELECTION,
        temperature_range_k: [300.0, 500.0],
    };
    let changed = Fixture::new(
        [900.0, 1500.0],
        [600.0, 1000.0],
        400.0,
        81000.0,
        [true, true],
    );
    assert!(IsobaricCrystallizationPair::try_new(&changed.card, query).is_err());
    let wrong_pressure = point(400.0, 200000.0);
    assert!(
        IsobaricCrystallizationPair::try_new(
            &fixture.card,
            IsobaricCrystallizationQuery {
                reference: &wrong_pressure,
                ..query
            }
        )
        .is_err()
    );
    assert!(
        IsobaricCrystallizationPair::try_new(
            &fixture.card,
            IsobaricCrystallizationQuery {
                temperature_range_k: [299.0, 500.0],
                ..query
            }
        )
        .is_err()
    );
    // A reference outside the kinetic interval remains valid when both sources cover it.
    let clipped = IsobaricCrystallizationPair::try_new(
        &fixture.card,
        IsobaricCrystallizationQuery {
            temperature_range_k: [300.0, 380.0],
            ..query
        },
    )
    .unwrap();
    close(clipped.latent_enthalpy_j_kg(350.0).unwrap(), 61250.0, 1e-10);
    assert!(clipped.parent().enthalpy_at_temperature(400.0).is_err());
}

fn thermal(law: &dyn ThermalLaw) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(1),
                segment_id: 1,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: 2.0,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [0.1, 0.0, 0.0],
                    volume_m3: 0.002,
                    surface_area_m2: 0.1,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        }],
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

#[test]
fn g1_native_crystallization_uses_exact_cp_recalescence_and_one_energy_account() {
    let fixture = Fixture::standard();
    let parent = fixture.caloric(&fixture.parent_key, 300.0);
    let crystal = fixture.caloric(&fixture.crystal_key, 300.0);
    let pair = fixture.pair(&parent, &crystal).unwrap();
    let rate = |_: f64| -> Result<f64, CrystallizationError> { Ok(0.7) };
    let rate_id =
        DomainHasher::new("org.frankensim.synthetic.constant-crystallization-rate").finalize();
    let law = pair.law(&rate, rate_id, 1.0, 1.0, 550.0).unwrap();
    let material = CrystallizationMaterial {
        material_id: MaterialId(1),
        amorphous: pair.parent(),
        law: &law,
        maximum_temperature_slope: pair.maximum_temperature_slope(),
    };
    let initial =
        CrystallizationSimulation::new(thermal(pair.parent()), &[material], || Ok(())).unwrap();
    let controls = CrystallizationControls {
        maximum_dt_s: 0.1,
        caloric: CrystallizationSolve {
            maximum_iterations: 128,
            specific_enthalpy_tolerance_j_kg: 1e-7,
            temperature_tolerance_k: 1e-10,
        },
        cell_energy_tolerance_j: 1e-6,
        total_energy_tolerance_j: 1e-4,
    };
    let (next, audit) = initial.advance_to(1.0, controls, || Ok(())).unwrap();
    let state = next.states()[&CellId(1)];
    // n=1 and constant K give X=1-exp(-Kt). Adiabatic enthalpy remains zero:
    // (1.5-0.5X)*(T²-400²) - 80000X = 0.
    let fraction = -(-0.7_f64).exp_m1();
    let expected_t = (160000.0 + 80000.0 * fraction / (1.5 - 0.5 * fraction)).sqrt();
    close(state.crystalline_mass_fraction(), fraction, 2e-13);
    close(state.temperature_k(), expected_t, 2e-9);
    close(state.specific_enthalpy_j_kg(), 0.0, 0.0);
    close(audit.thermal.external_heat_j, 0.0, 0.0);
    close(audit.thermal.stored_energy_change_j, 0.0, 0.0);
    close(audit.thermal.residual_j, 0.0, 0.0);
    close(audit.crystallized_mass_kg, 2.0 * fraction, 5e-13);
    assert!(state.temperature_k() > 430.0);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(
        initial.states()[&CellId(1)].crystalline_mass_fraction(),
        0.0
    );

    let failed = initial.advance_to(1.0, controls, || Err(DepositionError::Cancelled));
    assert!(failed.is_err());
    assert_eq!(initial.time_s(), 0.0);
    // The original unrebased Cp source cannot impersonate the registered phase chart.
    assert!(
        CrystallizationSimulation::new(
            thermal(&parent),
            &[CrystallizationMaterial {
                amorphous: &parent,
                ..material
            }],
            || Ok(())
        )
        .is_err()
    );
}

#[test]
fn g1_prepared_phase_integrals_cross_source_knots_and_preserve_clipped_endpoints() {
    let parent_key = PropertyKey::new("indexed_parent_cp", SPECIFIC_HEAT);
    let crystal_key = PropertyKey::new("indexed_crystal_cp", SPECIFIC_HEAT);
    let gap_key = PropertyKey::new("indexed_phase_gap", SPECIFIC_ENERGY);
    let temperature = QuantitySpec::dimensional(Temperature::DIMS);
    let pressure = QuantitySpec::dimensional(Pressure::DIMS);
    let provenance = || Provenance {
        source: "synthetic piecewise-Cp integration oracle".into(),
        license: "test-only".into(),
        artifact: None,
    };
    let mut claims = ClaimSet::new();
    for (key, capacities) in [
        (&parent_key, [900.0, 1050.0, 1200.0, 1500.0, 1800.0]),
        (&crystal_key, [600.0, 700.0, 800.0, 700.0, 600.0]),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: [300.0, 350.0, 400.0, 450.0, 500.0]
                        .into_iter()
                        .zip(capacities)
                        .collect(),
                    dims: SPECIFIC_HEAT,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", temperature, 300.0, 500.0)
                    .with_quantity("p", pressure, 100000.0, 100000.0),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    claims
        .insert_claim(PropertyClaim {
            key: gap_key.clone(),
            value: PropertyValue::Scalar {
                value: 80000.0,
                dims: SPECIFIC_ENERGY,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", temperature, 400.0, 400.0)
                .with_quantity("p", pressure, 100000.0, 100000.0),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    let fixture = Fixture {
        card: MaterialCard::assemble(
            MaterialStateId {
                chemistry: "one synthetic composition".into(),
                phase: "parent and crystal".into(),
                process: "piecewise-Cp numerical fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        parent_key,
        crystal_key,
        gap_key,
        reference: point(400.0, 100000.0),
    };
    let parent = fixture.caloric(&fixture.parent_key, 350.0);
    let crystal = fixture.caloric(&fixture.crystal_key, 425.0);
    let pair = IsobaricCrystallizationPair::try_new(
        &fixture.card,
        IsobaricCrystallizationQuery {
            parent: &parent,
            crystal: &crystal,
            enthalpy_gap: &fixture.gap_key,
            reference: &fixture.reference,
            selection: SELECTION,
            temperature_range_k: [325.0, 475.0],
        },
    )
    .unwrap();
    for branch in [pair.parent(), pair.crystal()] {
        assert!(branch.inverse.is_some());
        assert_eq!(branch.source().source_line().knots.len(), 5);
        assert_ne!(branch.source().datum_temperature_k(), 400.0);
        let span = branch.specific_enthalpy_span_j_kg();
        for (t, h) in [(325.0, span[0]), (475.0, span[1])] {
            assert_eq!(branch.enthalpy_at_temperature(t).unwrap(), h);
            assert_eq!(branch.temperature_k(h).unwrap(), t);
        }
        assert!(
            branch
                .enthalpy_at_temperature(325.0_f64.next_down())
                .is_err()
        );
        assert!(branch.enthalpy_at_temperature(475.0_f64.next_up()).is_err());
    }
    // The reference is 400 K. Queries at 337.5 and 462.5 K lie in different
    // source segments on BOTH sides, exercising prepared prefix+local work.
    // Clipped endpoints and their nearest interior values retain round trips.
    for t in [
        325.0,
        325.0_f64.next_up(),
        337.5,
        350.0,
        375.0,
        400.0,
        425.0,
        450.0,
        462.5,
        475.0_f64.next_down(),
        475.0,
    ] {
        // Cp_parent=3T below 400, then 1200+6(T-400).
        // Cp_crystal=2T below 400, then 800-2(T-400).
        let x = t - 400.0;
        let parent_h = if t <= 400.0 {
            1.5 * (t * t - 160000.0)
        } else {
            1200.0 * x + 3.0 * x * x
        };
        let crystal_h = -80000.0
            + if t <= 400.0 {
                t * t - 160000.0
            } else {
                800.0 * x - x * x
            };
        for (branch, expected) in [(pair.parent(), parent_h), (pair.crystal(), crystal_h)] {
            let h = branch.enthalpy_at_temperature(t).unwrap();
            close(h, expected, 3e-10);
            close(branch.temperature_k(h).unwrap(), t, 2e-12);
            close(branch.temperature_k(expected).unwrap(), t, 2e-12);
        }
        close(
            pair.latent_enthalpy_j_kg(t).unwrap(),
            parent_h - crystal_h,
            3e-10,
        );
    }
}

