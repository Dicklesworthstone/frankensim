#![cfg(feature = "enthalpy-coupling")]
//! G1/G3/G4: phase-dependent conductivity, spatial heat, and original energy.

use fs_blake3::ContentHash;
use fs_conduction::enthalpy::{ThermalLaw, ThermalStorageBasis};
use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ChtError, ConductivityTensor, EnergyConfig, FaceContactResistance, SolidMaterial,
    TENSOR_SDIRK2_GAMMA, TensorCaloricConductivity, TensorCaloricConductivityPoint,
    TensorCaloricEndpoint, TensorCaloricLaw, TensorCaloricStorage, TensorConductivityLaw,
    TensorNonlinearConfig, ThermalFace, ThermalSetup, Voxel, VoxelDomain,
    advance_caloric_tensor_conduction, advance_caloric_tensor_conduction_sdirk2,
    solve_nonlinear_tensor_conduction,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};

const MASS: f64 = 2.0;
const OPENING_H: f64 = 20.0;

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-7,
        "{actual:e} != {expected:e}"
    );
}

fn controls() -> TensorNonlinearConfig {
    TensorNonlinearConfig {
        linear: EnergyConfig {
            tolerance: 1e-10,
            ..EnergyConfig::default()
        },
        max_newton_iterations: 64,
        max_backtracks: 128,
    }
}

fn isotropic(value: f64) -> ConductivityTensor {
    ConductivityTensor::new([[value, 0.0, 0.0], [0.0, value, 0.0], [0.0, 0.0, value]]).unwrap()
}

fn melting(datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0xF1; 32]),
        [
            (-200.0, 200.0, 0.0),
            (0.0, 300.0, 0.0),
            (100.0, 300.0, 1.0),
            (500.0, 500.0, 1.0),
        ]
        .map(
            |(h, temperature_k, liquid_mass_fraction)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: datum + h,
                temperature_k,
                liquid_mass_fraction,
                bulk_density_kg_m3: MASS,
            },
        )
        .to_vec(),
    )
    .unwrap()
}

fn sensible(datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0xF2; 32]),
        SolidLiquidPhase::Solid,
        [(-200.0, 200.0), (400.0, 500.0)]
            .map(|(h, temperature_k)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: datum + h,
                temperature_k,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: MASS,
            })
            .to_vec(),
    )
    .unwrap()
}

/// Declared analytic test law k=1+alpha*(h-datum). The original caloric
/// identity and h coordinate are mandatory: Kelvin cannot select its phase.
struct PhaseConductivity<'a> {
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    datum: f64,
    alpha: f64,
    refuse_above: Option<f64>,
    cancel_above: Option<(f64, &'a CancelGate)>,
}

impl<'a> PhaseConductivity<'a> {
    fn new(curve: &'a EquilibriumEnthalpyPhaseCurve, datum: f64, alpha: f64) -> Self {
        Self {
            curve,
            datum,
            alpha,
            refuse_above: None,
            cancel_above: None,
        }
    }
}

impl TensorConductivityLaw for PhaseConductivity<'_> {
    fn conductivity_at_temperature(&self, _: f64) -> Result<ConductivityTensor, ChtError> {
        Err(ChtError::InvalidInput {
            field: "test.phase_kelvin",
            reason: "temperature alone does not select the original latent inventory".into(),
        })
    }

    fn mean_conductivity(&self, _: f64, _: f64) -> Result<ConductivityTensor, ChtError> {
        self.conductivity_at_temperature(300.0)
    }

    fn caloric_conductivity(
        &self,
        point: TensorCaloricConductivityPoint,
    ) -> Result<Option<TensorCaloricConductivity>, ChtError> {
        if point.storage_basis != ThermalStorageBasis::Enthalpy
            || point.caloric_identity != self.curve.caloric_identity()
        {
            return Err(ChtError::InvalidInput {
                field: "thermal.tensor_caloric.original_source_refusal",
                reason: "phase conductivity has a different original caloric owner".into(),
            });
        }
        let h = point.specific_energy_j_kg - self.datum;
        if self.refuse_above.is_some_and(|limit| h > limit) {
            return Err(ChtError::InvalidInput {
                field: "thermal.tensor_caloric.original_source_refusal",
                reason: "original phase-conductivity source refused a later query".into(),
            });
        }
        if let Some((limit, gate)) = self.cancel_above {
            if h > limit {
                gate.request();
                return Err(ChtError::Cancelled);
            }
        }
        near(
            point.temperature_k,
            self.curve
                .temperature_k(point.specific_energy_j_kg)
                .unwrap(),
        );
        TensorCaloricConductivity::new(
            isotropic(1.0 + self.alpha * h),
            [
                [self.alpha, 0.0, 0.0],
                [0.0, self.alpha, 0.0],
                [0.0, 0.0, self.alpha],
            ],
        )
        .map(Some)
    }
}

fn hot_face() -> (VoxelDomain, [SolidMaterial; 1], ThermalSetup) {
    let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("declared phase conductivity", 1.0)];
    let mut boundaries = [ThermalFace::Adiabatic; 6];
    boundaries[0] = ThermalFace::Temperature(310.0);
    (domain, solids, ThermalSetup::new(boundaries))
}

#[test]
fn latent_heating_changes_conductivity_at_constant_kelvin_and_preserves_the_energy_datum() {
    let (domain, solids, setup) = hot_face();
    let saved = setup.clone();
    for datum in [0.0, 1e6] {
        let curve = melting(datum);
        let k = PhaseConductivity::new(&curve, datum, 0.01);
        for dt in [0.25, 1.0] {
            let previous = [datum + OPENING_H];
            // A unit-area half cell gives q_in=20*k(h) W. Thus the exact
            // implicit residual has derivative M/dt-0.2 even though dT/dh=0.
            let expected_h = (OPENING_H + 10.0 * dt) / (1.0 - 0.1 * dt);
            let step = advance_caloric_tensor_conduction(
                &domain,
                &solids,
                &setup,
                &[&k],
                TensorCaloricStorage {
                    laws: &[&curve],
                    reference_mass_kg: &[MASS],
                    previous_specific_energy_j_kg: &previous,
                },
                &[],
                0.0,
                dt,
                &controls(),
                &CancelGate::new(),
            )
            .unwrap();
            assert_eq!(step.endpoint.solution.solution.temperature, vec![300.0]);
            near(step.endpoint.specific_energy_j_kg[0] - datum, expected_h);
            near(
                curve
                    .values_at_specific_enthalpy(step.endpoint.specific_energy_j_kg[0])
                    .unwrap()
                    .liquid_mass_fraction,
                expected_h / 100.0,
            );
            let energy = MASS * (expected_h - OPENING_H);
            near(step.record.stored_energy_change_j, energy);
            near(step.record.boundary_outflow_j, -energy);
            near(step.record.source_j, 0.0);
            near(step.record.closure_j, 0.0);
            near(
                step.endpoint.solution.solution.outward_heat_w[0][0],
                -20.0 * (1.0 + 0.01 * expected_h),
            );
            assert_eq!(previous, [datum + OPENING_H]);
        }
    }
    assert_eq!(setup, saved);

    // The continuous latent evolution is h(t)+100=120*exp(0.1*t).
    // Compare actual repeated native BE steps with that independent solution.
    let curve = melting(0.0);
    let k = PhaseConductivity::new(&curve, 0.0, 0.01);
    let exact = 120.0 * 0.2_f64.exp() - 100.0;
    let mut errors = Vec::new();
    for count in [4, 8, 16] {
        let dt = 2.0 / f64::from(count);
        let mut energy = [OPENING_H];
        let mut inward_j = 0.0;
        for index in 0..count {
            let step = advance_caloric_tensor_conduction(
                &domain,
                &solids,
                &setup,
                &[&k],
                TensorCaloricStorage {
                    laws: &[&curve],
                    reference_mass_kg: &[MASS],
                    previous_specific_energy_j_kg: &energy,
                },
                &[],
                f64::from(index) * dt,
                dt,
                &controls(),
                &CancelGate::new(),
            )
            .unwrap();
            inward_j -= step.record.boundary_outflow_j;
            energy[0] = step.endpoint.specific_energy_j_kg[0];
        }
        near(MASS * (energy[0] - OPENING_H), inward_j);
        errors.push((energy[0] - exact).abs());
    }
    assert!(errors[0] > 0.1);
    for pair in errors.windows(2) {
        assert!((1.9..2.2).contains(&(pair[0] / pair[1])));
    }
}

/// Two unit half cells plus a .25 m²K/W contact. The receiver's k=1+.02h
/// and the donor's k=2 give R=.5/(1.4+.01J)+.5 after transfer J joules.
/// The donor has M=2 and Cp=2, hence its temperature is 340-J/4 K.
fn transferred_j(interval: f64, carried_j: f64) -> f64 {
    let a = 0.005 + 0.0025 * interval;
    let b = 1.2 - 0.005 * carried_j - 0.05 * interval;
    let c = 1.2 * carried_j + 56.0 * interval;
    2.0 * c / (b + (b * b + 4.0 * a * c).sqrt())
}

fn check_pair(endpoint: &TensorCaloricEndpoint, transfer: f64) {
    let h = [OPENING_H + transfer / MASS, 80.0 - transfer / MASS];
    let k = 1.0 + 0.02 * h[0];
    let hot_t = 300.0 + h[1] / 2.0;
    let heat = (hot_t - 300.0) / (0.5 / k + 0.5);
    for (&actual, expected) in endpoint.specific_energy_j_kg.iter().zip(h) {
        near(actual, expected);
    }
    let field = &endpoint.solution.solution;
    assert_eq!(field.temperature[0], 300.0);
    near(field.temperature[1], hot_t);
    near(field.outward_heat_w[0][1], -heat);
    near(field.outward_heat_w[1][0], heat);
    near(field.face_temperature_k[0][1], 300.0 + 0.5 * heat / k);
    near(field.face_temperature_k[1][0], hot_t - 0.25 * heat);
    near(
        field.face_temperature_k[1][0] - field.face_temperature_k[0][1],
        0.25 * heat,
    );
    near(MASS * (h[0] + h[1] - OPENING_H - 80.0), 0.0);
}

#[test]
fn spatial_contact_and_both_sdirk_stages_use_the_actual_latent_conductivity() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let solids = [SolidMaterial::new("two original thermal inventories", 1.0)];
    let receiver = melting(0.0);
    let donor = sensible(0.0);
    let phase_k = PhaseConductivity::new(&receiver, 0.0, 0.02);
    let donor_k = isotropic(2.0);
    let caloric: [&dyn TensorCaloricLaw; 2] = [&receiver, &donor];
    let previous = [OPENING_H, 80.0];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.face_contacts.push(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: 0.25,
    });
    for dt in [0.5, 3.0] {
        let storage = TensorCaloricStorage {
            laws: &caloric,
            reference_mass_kg: &[MASS; 2],
            previous_specific_energy_j_kg: &previous,
        };
        let backward = advance_caloric_tensor_conduction(
            &domain,
            &solids,
            &setup,
            &[&phase_k, &donor_k],
            storage,
            &[],
            0.0,
            dt,
            &controls(),
            &CancelGate::new(),
        )
        .unwrap();
        check_pair(&backward.endpoint, transferred_j(dt, 0.0));
        near(backward.record.stored_energy_change_j, 0.0);
        near(backward.record.boundary_outflow_j, 0.0);
        near(backward.record.closure_j, 0.0);

        let a = TENSOR_SDIRK2_GAMMA * dt;
        let b = (1.0 - TENSOR_SDIRK2_GAMMA) * dt;
        let first_transfer = transferred_j(a, 0.0);
        let staged = advance_caloric_tensor_conduction_sdirk2(
            &domain,
            &solids,
            [&setup, &setup],
            &[&phase_k, &donor_k],
            storage,
            [&[], &[]],
            0.0,
            dt,
            &controls(),
            &CancelGate::new(),
        )
        .unwrap();
        check_pair(&staged.stages[0], first_transfer);
        check_pair(&staged.stages[1], transferred_j(a, b * first_transfer / a));
        near(staged.record.stored_energy_change_j, 0.0);
        near(staged.record.boundary_outflow_j, 0.0);
        near(staged.record.source_j, 0.0);
        near(staged.record.closure_j, 0.0);
    }
    assert_eq!(previous, [OPENING_H, 80.0]);
}

#[test]
fn phase_sources_refuse_kelvin_only_mismatched_and_late_failed_queries_atomically() {
    let (domain, solids, setup) = hot_face();
    let saved = setup.clone();
    let curve = melting(0.0);
    let k = PhaseConductivity::new(&curve, 0.0, 0.01);
    assert!(matches!(
        solve_nonlinear_tensor_conduction(
            &domain,
            &solids,
            &setup,
            &[&k],
            &[300.0],
            &controls(),
            &CancelGate::new(),
        ),
        Err(ChtError::InvalidInput {
            field: "test.phase_kelvin",
            ..
        })
    ));

    let previous = [OPENING_H];
    let wrong_curve = melting(1.0);
    let gate = CancelGate::new();
    let solve = |law: &dyn TensorCaloricLaw, conductivity: &dyn TensorConductivityLaw| {
        advance_caloric_tensor_conduction(
            &domain,
            &solids,
            &setup,
            &[conductivity],
            TensorCaloricStorage {
                laws: &[law],
                reference_mass_kg: &[MASS],
                previous_specific_energy_j_kg: &previous,
            },
            &[],
            0.0,
            1.0,
            &controls(),
            &gate,
        )
    };
    assert!(matches!(
        solve(&wrong_curve, &k),
        Err(ChtError::InvalidInput {
            field: "thermal.tensor_caloric.original_source_refusal",
            ..
        })
    ));
    let mut refused = PhaseConductivity::new(&curve, 0.0, 0.01);
    refused.refuse_above = Some(OPENING_H);
    assert!(matches!(
        solve(&curve, &refused),
        Err(ChtError::InvalidInput {
            field: "thermal.tensor_caloric.original_source_refusal",
            ..
        })
    ));
    let mut cancelled = PhaseConductivity::new(&curve, 0.0, 0.01);
    cancelled.cancel_above = Some((OPENING_H, &gate));
    assert!(matches!(
        solve(&curve, &cancelled),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(previous, [OPENING_H]);
    assert_eq!(setup, saved);
}

#[test]
fn original_source_adapter_crosses_a_latent_knot_without_losing_its_restricted_chart() {
    use fs_evidence::ValidityDomain;
    use fs_lbm::conjugate::IsotropicPhaseConductivity;
    use fs_matdb::{
        ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
        UncertaintyModel,
    };
    use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS as KD;
    use fs_material::thermal::conduction::phase::PreparedPhaseConductivity;
    use fs_material::thermal::source::ThermalMaterialPoint;
    use fs_qty::{Dims, QuantitySpec, Temperature};

    let mut claims = ClaimSet::new();
    let key = PropertyKey::new("k", KD);
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "f".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, 1.0), (0.25, 1.5), (1.0, 2.25)],
                dims: KD,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    200.0,
                    500.0,
                )
                .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured latent conductivity fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured material".into(),
            phase: "solid-liquid".into(),
            process: "declared phase conductivity".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve =
        EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(), melting(0.0).knots().to_vec())
            .unwrap();
    let context = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap()
        .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.2)
        .unwrap();
    let source = PreparedPhaseConductivity::new(
        &curve,
        ThermalMaterialPoint {
            material: &card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: Some("f"),
        },
        &key,
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        [10.0, 100.0],
        || Ok(()),
    )
    .unwrap();
    let gate = CancelGate::new();
    let law = IsotropicPhaseConductivity::new(&source, &gate).unwrap();
    assert!(std::ptr::eq(law.source(), &source));
    assert_eq!(law.caloric_identity(), curve.caloric_identity());
    assert_eq!(law.specific_energy_bounds_j_kg(), [10.0, 100.0]);
    assert_eq!(
        source.knots().iter().map(|k| k[0]).collect::<Vec<_>>(),
        vec![10.0, 25.0, 100.0]
    );
    // The prepared restriction ends on the latent side of the original melt
    // boundary. Its tangent must not silently switch to the sensible branch.
    assert_eq!(law.temperature_slope_kg_k_j(100.0).unwrap(), 0.0);
    assert_eq!(curve.temperature_slope_kg_k_j(100.0).unwrap(), 0.5);
    assert!(law.temperature_k(100.0_f64.next_up()).is_err());
    assert!(curve.temperature_k(100.0_f64.next_up()).is_ok());
    assert!(law.conductivity_at_temperature(300.0).is_err());
    assert!(law.mean_conductivity(300.0, 310.0).is_err());

    let point = TensorCaloricConductivityPoint {
        specific_energy_j_kg: 100.0,
        temperature_k: 300.0,
        storage_basis: ThermalStorageBasis::Enthalpy,
        caloric_identity: curve.caloric_identity(),
    };
    let endpoint = law.caloric_conductivity(point).unwrap().unwrap();
    assert_eq!(endpoint.conductivity(), isotropic(2.25));
    near(endpoint.derivative_kg_w_m_k_j()[0][0], 0.01);
    for invalid in [
        TensorCaloricConductivityPoint {
            storage_basis: ThermalStorageBasis::InternalEnergy,
            ..point
        },
        TensorCaloricConductivityPoint {
            caloric_identity: melting(0.0).caloric_identity(),
            ..point
        },
        TensorCaloricConductivityPoint {
            temperature_k: 301.0,
            ..point
        },
    ] {
        assert!(matches!(
            law.caloric_conductivity(invalid),
            Err(ChtError::InvalidInput {
                field: "thermal.tensor_caloric.original_source_refusal",
                ..
            })
        ));
    }

    let (domain, solids, setup) = hot_face();
    let previous = [OPENING_H];
    let step = advance_caloric_tensor_conduction(
        &domain,
        &solids,
        &setup,
        &[&law],
        TensorCaloricStorage {
            laws: &[&law],
            reference_mass_kg: &[MASS],
            previous_specific_energy_j_kg: &previous,
        },
        &[],
        0.0,
        1.0,
        &controls(),
        &gate,
    )
    .unwrap();
    // Crossing h=25 changes dk/dh from .02 to .01. On the final segment,
    // 2*(h-20)=20*(1.25+.01*h), while T stays exactly 300 K.
    let expected_h = 65.0 / 1.8;
    near(step.endpoint.specific_energy_j_kg[0], expected_h);
    assert_eq!(step.endpoint.solution.solution.temperature, vec![300.0]);
    near(
        step.record.stored_energy_change_j,
        MASS * (expected_h - OPENING_H),
    );
    near(
        step.record.boundary_outflow_j,
        -MASS * (expected_h - OPENING_H),
    );
    near(step.record.closure_j, 0.0);
    assert_eq!(previous, [OPENING_H]);
}
