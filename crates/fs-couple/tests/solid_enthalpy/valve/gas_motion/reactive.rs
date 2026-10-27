//! G1/G3/G4: source chemistry drives the original mesh mode, with an independent
//! constant-cv pressure/work oracle and transactional refusal controls.
//! All coefficients are synthetic; these are not empirical combustion checks.
use super::*;
use fs_couple::cavity_phs::modal::reactive::ReactiveModalGasError;
use fs_couple::cavity_phs::piston::reactive::{
    GasReaction, GasReactionError, GasReactionSolve,
};
use fs_couple::cavity_phs::piston::reactive::material_rate::{
    ACTIVITY_RATE_SCALE_DIMS, ACTIVITY_RATE_SCALE_PROPERTY, MaterialRateScale,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, UncertaintyModel,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_thermochem::equilibrium::{
    DeclaredStoichiometricMassActionV1, IdealGasReactionEquilibriumV1,
    MassActionKineticOrderConventionV1, ReactionProgressRateQuantityV1,
    ReactionProgressRateScaleV1, ThermodynamicReverseRateClosureV1,
};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ReactionId, SpeciesId,
    StoichiometricMatrix, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, verify_conservation,
};

const RELEASE: f64 = 4000.0;
const TEMPERATURE: f64 = 400.0;

// Balanced isomers have identical constant cv and molar mass. Only B's
// formation reference differs. The source-rate card, not a reaction or motion
// branch, supplies the varied forward scale in mol/(m3 s).
fn reacting_gas(scale: f64, fractions: [f64; 3]) -> (GasReaction, ClosedGasParcel) {
    let mut components = Vec::new();
    let mut transport = Vec::new();
    for (name, formation) in [("A", 0.0), ("B", -RELEASE), ("I", 0.0)] {
        let (caloric, source) = gas_fixture::component(name, 0.028, formation);
        components.push(caloric);
        transport.push(source);
    }
    let catalog = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            components,
            Composition::new(CompositionBasis::MoleFraction, vec![0.25, 0.25, 0.5])
                .unwrap(),
        )
        .unwrap(),
        transport,
        None,
    )
    .unwrap();
    let species = vec![SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()];
    let element = ElementalMatrix::new(
        vec![ElementId::new("X").unwrap()],
        species.clone(),
        vec![vec![1, 1]],
    )
    .unwrap();
    let id = ReactionId::new("elastic-cavity-isomers").unwrap();
    let matrix = StoichiometricMatrix::new(
        species.clone(),
        vec![id.clone()],
        vec![vec![-1], vec![1]],
    )
    .unwrap();
    let charge = ChargeVector::new(species, vec![0, 0]).unwrap();
    let certificate = verify_conservation(&element, &matrix, &charge).unwrap();
    let equilibrium = IdealGasReactionEquilibriumV1::new(
        matrix,
        certificate,
        id,
        catalog.thermodynamics().components()[..2].to_vec(),
    )
    .unwrap();
    let kinetics = DeclaredStoichiometricMassActionV1::new(
        ThermodynamicReverseRateClosureV1::new(
            equilibrium,
            ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(1.0))
                .unwrap(),
        ),
        MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients,
    );
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(ACTIVITY_RATE_SCALE_PROPERTY, ACTIVITY_RATE_SCALE_DIMS),
            value: PropertyValue::Scalar {
                value: scale,
                dims: ACTIVITY_RATE_SCALE_DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic source-defined elastic-cavity kinetics".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let rate = MaterialRateScale::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic-isomers".into(),
                phase: "gas".into(),
                process: "authored".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        Pressure::new(EXTERIOR),
    )
    .unwrap();
    let reaction = GasReaction::new(catalog, kinetics)
        .unwrap()
        .with_material_rate(rate)
        .unwrap();
    let model = reaction
        .catalog()
        .with_composition(
            Composition::new(CompositionBasis::MoleFraction, fractions.to_vec()).unwrap(),
        )
        .unwrap();
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            TEMPERATURE,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), EXTERIOR)
        .unwrap();
    (reaction, ClosedGasParcel::new(model, point, VOLUME).unwrap())
}

fn chemistry() -> GasReactionSolve {
    GasReactionSolve {
        gas: GasEnergySolve {
            energy_tolerance_j_kg: 1e-7,
            ..controls().gas
        },
        maximum_consumed_fraction: 0.1,
        species_tolerance_mol: 1e-13,
        entropy_tolerance_j_k: 1e-11,
    }
}

struct Reference {
    volume_change: f64,
    velocity: f64,
    pressure: f64,
    temperature: f64,
}

// Direct integrals on the existing six-tet cube (rho=6, nu=0, node7.x free):
// M=rho/10, K=2E/3. Do not read the reduced mass, frequency, force map,
// gas work or pressure from the implementation being checked.
// All three species have cv=2.5R and equal molar mass; n is conserved.
// Chemistry changes the formation offset by -RELEASE*xi at conserved total U,
// hence T*=T0+RELEASE*xi/(n*cv). The native energy-midpoint gas equation is
// Pmid=nRT*/(V0+gamma*dV/2), not Pmid=(gamma-1)*TOTAL_U/V.
// Starting at rest, exact held-force motion gives dV=C*(Pmid-Pext).
fn reference(young: f64, extent: f64) -> Reference {
    let dt = 1.0 / RATE as f64;
    let mass = 6.0 / 10.0;
    let stiffness = 2.0 * young / 3.0;
    let omega = (stiffness / mass).sqrt();
    let compliance = AREA * AREA * (1.0 - (omega * dt).cos()) / stiffness;
    let amount = EXTERIOR * VOLUME / (R * TEMPERATURE);
    let temperature = TEMPERATURE + RELEASE * extent / (amount * 2.5 * R);
    let drive = amount * R * temperature - EXTERIOR * VOLUME;
    let b = VOLUME + 0.7 * compliance * EXTERIOR;
    // Stable root of .7*dV^2 + b*dV - C*drive = 0, for either reaction sign.
    let volume_change =
        2.0 * compliance * drive / (b + (b * b + 2.8 * compliance * drive).sqrt());
    let pressure = amount * R * temperature / (VOLUME + 0.7 * volume_change);
    let velocity = AREA * (pressure - EXTERIOR) * (omega * dt).sin() / (mass * omega);
    Reference {
        volume_change,
        velocity,
        pressure,
        temperature: temperature - pressure * volume_change / (amount * 2.5 * R),
    }
}

#[test]
fn source_defined_chemistry_drives_the_actual_mode_and_matches_independent_work() {
    with_cx(|cx| {
        let mut velocities = Vec::new();
        for (scale, young) in [(100.0, 1000.0), (300.0, 1000.0), (100.0, 2000.0)] {
            let (reaction, gas) = reacting_gas(scale, [1.0, 0.0, 0.0]);
            let initial_u = gas.state().specific_internal_energy_j_kg();
            let initial_mass = gas.mass_kg();
            let solid = motion(young, 0.0, 0.0, AREA, RATE, cx);
            let identity = solid.valve().identity();
            let mut run = GasValveMotion::new(gas, solid, cx).unwrap();
            let report = run
                .step_reactive(&reaction, 0.0, EXTERIOR, chemistry(), controls(), cx)
                .unwrap();
            let extent = scale * VOLUME / RATE as f64; // initial activity of A=1
            let expected = reference(young, extent);
            let physical = run.motion().nodal_motion(cx).unwrap();
            near(report.reaction.extent_mol, extent, 1e-13);
            assert!(report.reaction.rate_material.is_some());
            near(
                physical.displacement_m[7][0],
                expected.volume_change / AREA,
                1e-11,
            );
            near(physical.velocity_m_s[7][0], expected.velocity, 1e-10);
            near(
                report.mechanical.midpoint.state().pressure,
                expected.pressure,
                1e-6,
            );
            near(run.gas().state().state().temperature, expected.temperature, 1e-6);
            near(
                report.gas_energy_change_j,
                -expected.pressure * expected.volume_change,
                controls().energy_tolerance_j,
            );
            let measured = initial_mass
                * (run.gas().state().specific_internal_energy_j_kg() - initial_u)
                + energy(run.motion())
                + EXTERIOR * AREA * physical.displacement_m[7][0];
            assert!(measured.abs() <= controls().energy_tolerance_j);
            assert!(report.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
            assert_eq!(run.gas().mass_kg(), initial_mass);
            assert_eq!(run.motion().valve().identity(), identity);
            assert_eq!(run.motion().accepted_steps(), 1);
            assert_eq!(report.mechanical.motion.time_s, 1.0 / RATE as f64);
            velocities.push(physical.velocity_m_s[7][0]);
        }
        assert!(velocities[0] > 0.0);
        assert!(velocities[1] > 2.9 * velocities[0]);
        assert!(velocities[2] < velocities[0] - 1e-10);
    });
}

#[test]
fn reverse_chemistry_cools_and_closes_the_wall_without_a_negative_heat_preset() {
    with_cx(|cx| {
        let scale = 100.0;
        let (reaction, gas) = reacting_gas(scale, [0.0, 1.0, 0.0]);
        let mut run =
            GasValveMotion::new(gas, motion(1000.0, 0.0, 0.0, AREA, RATE, cx), cx)
                .unwrap();
        let report = run
            .step_reactive(&reaction, 0.0, EXTERIOR, chemistry(), controls(), cx)
            .unwrap();
        // Equal cp/entropy references give K=exp(RELEASE/(R*T)); initial B=1.
        let extent = -scale * VOLUME / RATE as f64 * (-RELEASE / (R * TEMPERATURE)).exp();
        let expected = reference(1000.0, extent);
        let physical = run.motion().nodal_motion(cx).unwrap();
        near(report.reaction.extent_mol, extent, 1e-13);
        near(physical.velocity_m_s[7][0], expected.velocity, 1e-10);
        near(run.gas().state().state().temperature, expected.temperature, 1e-6);
        assert!(report.reaction.extent_mol < 0.0);
        assert!(physical.displacement_m[7][0] < 0.0);
        assert!(physical.velocity_m_s[7][0] < 0.0);
        assert!(run.gas().state().state().temperature < TEMPERATURE);
        assert!(report.gas_energy_change_j > 0.0); // exterior compresses the cooled gas
        assert!(report.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
    });
}

#[test]
fn absent_reagents_keep_the_native_path_and_continuing_reactions_retain_the_anchor() {
    with_cx(|cx| {
        let (reaction, gas) = reacting_gas(100.0, [0.0, 0.0, 1.0]);
        let mut reactive =
            GasValveMotion::new(gas, motion(1000.0, 0.02, 0.1, AREA, RATE, cx), cx)
                .unwrap();
        let mut native = reactive.clone();
        for tick in 1..=4 {
            let report = reactive
                .step_reactive(&reaction, 0.0, EXTERIOR, chemistry(), controls(), cx)
                .unwrap();
            native.step(0.0, EXTERIOR, controls(), cx).unwrap();
            assert_eq!(report.reaction.extent_mol, 0.0);
            assert_eq!(reactive.gas(), native.gas());
            assert_eq!(reactive.motion().model().states(), native.motion().model().states());
            assert_eq!(reactive.motion().accepted_steps(), tick);
        }

        let (reaction, gas) = reacting_gas(100.0, [1.0, 0.0, 0.0]);
        let initial_u = gas.state().specific_internal_energy_j_kg();
        let mass = gas.mass_kg();
        let mut run =
            GasValveMotion::new(gas, motion(1000.0, 0.0, 0.1, AREA, RATE, cx), cx)
                .unwrap();
        let mut physical_loss = 0.0;
        for tick in 1..=8 {
            let report = run
                .step_reactive(&reaction, 0.0, EXTERIOR, chemistry(), controls(), cx)
                .unwrap();
            let physical = run.motion().nodal_motion(cx).unwrap();
            near(
                run.gas().volume_m3() - VOLUME,
                AREA * physical.displacement_m[7][0],
                controls().volume_tolerance_m3,
            );
            physical_loss += report.mechanical.motion.modal.viscous_dissipation_j;
            let measured = mass * (run.gas().state().specific_internal_energy_j_kg() - initial_u)
                + energy(run.motion())
                + physical_loss
                + EXTERIOR * AREA * physical.displacement_m[7][0];
            assert!(measured.abs() <= tick as f64 * controls().energy_tolerance_j);
            assert!(report.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
            assert_eq!(run.motion().accepted_steps(), tick);
            assert!(report.reaction.extent_mol > 0.0);
        }
        assert!(physical_loss > 0.0);
    });
}

#[test]
fn chemical_and_late_mechanical_refusals_leave_all_histories_retryable() {
    with_cx(|cx| {
        let (reaction, gas) = reacting_gas(100.0, [1.0, 0.0, 0.0]);
        let mut run =
            GasValveMotion::new(gas, motion(1000.0, 0.0, 0.0, AREA, RATE, cx), cx)
                .unwrap();
        let original = run.clone();
        let unchanged = |run: &GasValveMotion| {
            assert_eq!(run.gas(), original.gas());
            assert_eq!(run.motion().model().states(), original.motion().model().states());
            assert_eq!(run.motion().accepted_steps(), original.motion().accepted_steps());
            assert_eq!(run.motion().valve().identity(), original.motion().valve().identity());
        };
        let mut limited_chemistry = chemistry();
        limited_chemistry.maximum_consumed_fraction = 1e-8;
        assert!(matches!(
            run.step_reactive(&reaction, 0.0, EXTERIOR, limited_chemistry, controls(), cx),
            Err(ReactiveModalGasError::Reaction(GasReactionError::ConsumptionLimit { .. }))
        ));
        unchanged(&run);
        let mut limited_motion = controls();
        limited_motion.maximum_trials = 1;
        assert!(matches!(
            run.step_reactive(&reaction, 0.0, EXTERIOR, chemistry(), limited_motion, cx),
            Err(ReactiveModalGasError::Motion(ModalGasError::TrialBudget))
        ));
        unchanged(&run);
        let preview = run
            .reactive_step_candidate(&reaction, 0.0, EXTERIOR, chemistry(), controls(), cx)
            .unwrap();
        unchanged(&run);
        let report = run
            .step_reactive(&reaction, 0.0, EXTERIOR, chemistry(), controls(), cx)
            .unwrap();
        assert_eq!(run.gas(), preview.motion().gas());
        assert_eq!(run.motion().model().states(), preview.motion().motion().model().states());
        assert_eq!(report.reaction, preview.step().reaction);
        assert_eq!(run.motion().accepted_steps(), 1);
    });
}
