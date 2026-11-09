//! The original molecular source must survive both layer and volume feedback.
use super::*;
use crate::cavity_modal::interior::heterogeneous::porous::AcousticVolumeRegionState;
use crate::porous_heating::feedback::{IsobaricRelaxingPorousGas, PorousRelaxingState};
use fs_duct::porous::PorousLayer;
use fs_material::MaterialError;
use fs_material::porous::PorousFluidResponse;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation, VibrationalMode,
};
use fs_thermochem::SpeciesId;

fn molecular_source() -> (RelaxationAbsorption, RelaxingMedium) {
    let species = SpeciesId::new("synthetic").unwrap();
    let law = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes: (0..2)
            .map(|i| RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: [1300.0, 1900.0][i],
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 101325.0 * [0.0002, 0.003][i],
                    t_ref_k: 300.0,
                    exponent: 1.0,
                    t_min_k: 250.0,
                    t_max_k: 550.0,
                },
            })
            .collect(),
    }])
    .unwrap();
    let reference = RelaxingMedium::new(
        GasState::try_new(&SPEC, 300.0, 101325.0).unwrap(),
        vec![(species, 1.0)],
    )
    .unwrap();
    (law, reference)
}

#[test]
fn accepted_temperatures_preserve_original_molecular_free_and_porous_volume_response() {
    let (mesh, positions) = fixture();
    let geometry = geometry(&mesh, &positions);
    let caloric = caloric();
    let porous = laws();
    let (relaxation, reference) = molecular_source();
    let source = IsobaricRelaxingPorousGas::try_new(
        IsobaricPorousGas::try_new(&SPEC, 101325.0, [280.0, 500.0]).unwrap(),
        &relaxation,
        &reference,
    )
    .unwrap();
    let thermal = initial(events([0.0, 40.0 * CP], [0.12, 0.1], false), &caloric);
    let mut destinations = mapping(&caloric, Some(&porous));
    destinations[0].porous = None;
    let problem = problem(false);
    let heating = TemperatureInteriorHeating::try_new(
        &geometry,
        &thermal,
        &destinations,
        &[],
        source,
        problem.clone(),
        2,
        || Ok(()),
    )
    .unwrap();
    let initial_reservoirs: Vec<_> = heating
        .fluid_states()
        .iter()
        .map(|s| s.source().clone())
        .collect();
    let result = heating
        .advance_to(&thermal, 0.0625, config(0.03125), || Ok(()))
        .unwrap();
    assert_eq!(result.work.fluid_evaluations, 4);
    assert_eq!(result.work.porous_evaluations, 2);
    assert_eq!(result.work.field_solves, 2);

    for (owner, native_thermal) in [
        (&heating, &thermal),
        (&result.heating, &result.step.thermal),
    ] {
        assert!(core::ptr::eq(owner.source().relaxation_law(), &relaxation));
        assert!(core::ptr::eq(owner.source().reference_medium(), &reference));
        // Independent orchestration reevaluates the ORIGINAL reservoir law at
        // each native temperature, then uses the original mixed-region adapter.
        let states = IDS.map(|cell| {
            let gas =
                GasState::try_new(&SPEC, native_thermal.temperature_k(cell).unwrap(), 101325.0)
                    .unwrap();
            relaxation
                .evaluate(&RelaxingMedium::new(gas, reference.components().to_vec()).unwrap())
                .unwrap()
        });
        let sources = [
            AcousticVolumeRegion::MolecularGas(&states[0]),
            AcousticVolumeRegion::RelaxingPorous {
                law: &porous[1],
                relaxation: &states[1],
            },
        ];
        let original = PorousRegionsInterior::new(
            &geometry,
            FREQUENCY,
            &sources,
            &[0, 1],
            &[],
            &[],
            problem.controls,
            &mut || false,
        )
        .unwrap();
        let field = original
            .operator()
            .solve(&problem.inward_velocity, &mut || false)
            .unwrap();
        assert_eq!(owner.field().pressure_pa, field.pressure_pa);
        assert_eq!(owner.field().element_power, field.element_power);
        for (i, state) in states.iter().enumerate() {
            let cached = &owner.fluid_states()[i];
            assert_eq!(cached.source(), state);
            let cloned = cached.clone();
            assert!(core::ptr::eq(cached.source(), cloned.source()));
            assert_eq!(
                owner.material_responses()[i],
                original.region_state(i).unwrap()
            );
        }
        let omega = core::f64::consts::TAU * FREQUENCY;
        let free = PorousFluidResponse::relaxing_fluid(&states[0], omega).unwrap();
        assert_eq!(
            owner.material_responses()[0],
            AcousticVolumeRegionState::MolecularGas(free)
        );
        assert_eq!(free.flow_resistivity_pa_s_m2(), 0.0);
        assert!(free.dynamic_bulk_modulus_pa().im > 0.0);
        let joint = porous[1].effective_relaxing(&states[1], omega).unwrap();
        assert_eq!(
            owner.material_responses()[1],
            AcousticVolumeRegionState::Porous(joint)
        );
        assert_ne!(
            joint.dynamic_bulk_modulus_pa(),
            porous[1]
                .effective_fluid(states[1].gas(), omega)
                .unwrap()
                .dynamic_bulk_modulus_pa()
        );
        // The pre-existing layer hook retains its same original joint response.
        let layer = source
            .prepare_layer(
                &porous[1],
                &owner.fluid_states()[1],
                omega,
                0.02,
                1e-4,
                states[1].gas().characteristic_impedance / 1e-4,
                0.0,
            )
            .unwrap();
        assert_eq!(*layer.material_response(), joint);
        assert!(core::ptr::eq(
            layer.relaxation_state().unwrap(),
            owner.fluid_states()[1].source()
        ));
        assert_eq!(field.element_power[0].density_dissipation_w, 0.0);
        assert!(field.element_power[0].compressional_dissipation_w > 0.0);
        assert!(field.element_power[1].density_dissipation_w > 0.0);
        assert!(field.power_balance_residual_w().abs() < 1e-9);
    }
    for (i, initial) in initial_reservoirs.iter().enumerate() {
        assert_eq!(heating.fluid_states()[i].source(), initial);
        assert_ne!(
            result.heating.fluid_states()[i].source().processes(),
            initial.processes()
        );
        assert_eq!(
            result.step.thermal.cells()[&IDS[i]].mass_kg,
            thermal.cells()[&IDS[i]].mass_kg
        );
    }
    let heat = result.step.acoustic;
    assert_eq!(heat.wall_resistance_dissipation_j, 0.0);
    assert_eq!(heat.wall_relaxation_dissipation_j, 0.0);
    close(
        heat.deposited_j,
        heat.density_dissipation_j + heat.compressional_dissipation_j,
        2e-12,
    );
    close(
        result.step.heat_report.external_heat_j,
        heat.deposited_j + 0.22 * 0.0625,
        2e-11,
    );
    close(result.step.heat_report.residual_j, 0.0, 2e-11);
}

/// A source authored for the old stack API, with additional constitutive state.
/// Its equilibrium ThermoviscousFluid implementation cannot authorize a volume.
#[derive(Clone, Copy)]
struct StackOnly<'a>(IsobaricRelaxingPorousGas<'a>);
impl PorousTemperatureFluid for StackOnly<'_> {
    type State = PorousRelaxingState;
    fn evaluate(&self, temperature: f64) -> Result<Self::State, PorousHeatingError> {
        self.0.evaluate(temperature)
    }
    fn temperature_k(&self, state: &Self::State) -> f64 {
        self.0.temperature_k(state)
    }
    #[allow(clippy::too_many_arguments)]
    fn prepare_layer<'a>(
        &self,
        material: &'a PorousFluidLaw,
        state: &Self::State,
        omega: f64,
        thickness_m: f64,
        area_m2: f64,
        reference_impedance: f64,
        tangential_wavenumber_rad_m: f64,
    ) -> Result<PorousLayer<'a>, PorousHeatingError> {
        self.0.prepare_layer(
            material,
            state,
            omega,
            thickness_m,
            area_m2,
            reference_impedance,
            tangential_wavenumber_rad_m,
        )
    }
}

#[test]
fn layer_only_dispersive_source_refuses_undeclared_volume_instead_of_losing_reservoirs() {
    let (mesh, positions) = fixture();
    let geometry = geometry(&mesh, &positions);
    let caloric = caloric();
    let porous = laws();
    let (relaxation, reference) = molecular_source();
    let source = StackOnly(
        IsobaricRelaxingPorousGas::try_new(
            IsobaricPorousGas::try_new(&SPEC, 101325.0, [280.0, 500.0]).unwrap(),
            &relaxation,
            &reference,
        )
        .unwrap(),
    );
    let state = source.evaluate(300.0).unwrap();
    let omega = core::f64::consts::TAU * FREQUENCY;
    assert!(
        source
            .prepare_layer(&porous[0], &state, omega, 0.02, 1e-4, 4e6, 0.0)
            .is_ok()
    );
    for material in [None, Some(&porous[0])] {
        assert!(matches!(
            source.volume_response(&state, material, omega),
            Err(MaterialError::Parameters { .. })
        ));
    }
    let thermal = initial(events([0.0; 2], [0.0; 2], false), &caloric);
    assert!(matches!(
        TemperatureInteriorHeating::try_new(
            &geometry,
            &thermal,
            &mapping(&caloric, Some(&porous)),
            &[],
            source,
            problem(false),
            2,
            || Ok(()),
        ),
        Err(InteriorFeedbackError::Material {
            element: 0,
            source: MaterialError::Parameters { .. }
        })
    ));
}
