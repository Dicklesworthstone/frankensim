//! Collective spherical geometry on the original shared gas donor and calorics.
use super::*;
use crate::evaporation::gas::droplet::{
    advance_droplet_with_gas_and_heat_at_sherwood, advance_droplet_with_gas_at_sherwood,
};
use crate::evaporation::gas::tests::{Phase, binding, catalog, controls, gas, liquid, temperature};
use fs_thermal::evaporation::{
    EvaporationProperties, InterfaceKinetics, KineticEvaporation, MOLAR_GAS_CONSTANT,
};
use fs_thermochem::SpeciesId;

const D: f64 = 2e-5;

fn row<'a>(
    id: u64,
    binding: &'a GasLiquidBinding,
    mass: f64,
    t: f64,
    d: f64,
    sh: f64,
) -> GasDropletSurface<'a> {
    GasDropletSurface {
        id,
        law: &Phase,
        binding,
        droplet: SphericalDropletState {
            mass_kg: mass,
            specific_enthalpy_j_kg: liquid(t).specific_enthalpy_j_kg,
        },
        temperature_from_enthalpy: &temperature,
        binary_diffusivity_m2_s: d,
        sherwood_number: sh,
    }
}
fn solve() -> GasLiquidChamberSolve {
    GasLiquidChamberSolve {
        surface: controls(),
        maximum_surfaces: 8,
    }
}
fn heat_controls() -> GasFilmHeatSolve {
    GasFilmHeatSolve {
        pure_liquid_enthalpy_range_j_kg: [
            liquid(200.0).specific_enthalpy_j_kg,
            liquid(900.0).specific_enthalpy_j_kg,
        ],
        temperature_tolerance_k: 1e-8,
        maximum_iterations: 128,
    }
}
fn outer() -> GasLiquidChamberHeatSolve {
    GasLiquidChamberHeatSolve {
        maximum_iterations: 128,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance {tolerance:e}"
    );
}
fn vapor() -> SpeciesId {
    SpeciesId::new("B").unwrap()
}
fn carrier() -> SpeciesId {
    SpeciesId::new("A").unwrap()
}

#[test]
fn different_sizes_diffusivities_and_sherwood_numbers_share_one_original_gas() {
    let source = catalog(-1e7);
    let binding = binding(-1e7);
    let original = gas(&source, 320.0, 1600.0, 0.01);
    let saved = original.clone();
    let inputs = [
        row(9, &binding, 8e-6, 280.0, D / 2.0, 2.0),
        row(2, &binding, 1e-6, 360.0, D, 5.0),
        row(5, &binding, 3e-6, 300.0, D * 1.5, 8.0),
    ];
    let dt = 0.02;
    let candidate =
        advance_droplets_with_gas(&inputs, &original, &source, dt, solve(), || false).unwrap();
    let reversed = advance_droplets_with_gas(
        &[inputs[2], inputs[1], inputs[0]],
        &original,
        &source,
        dt,
        solve(),
        || false,
    )
    .unwrap();
    assert_eq!(candidate, reversed);
    assert_eq!(
        candidate
            .droplets()
            .iter()
            .map(|x| x.id)
            .collect::<Vec<_>>(),
        [2, 5, 9]
    );

    let original_surfaces: Vec<_> = inputs
        .iter()
        .map(|input| {
            let t = temperature(input.droplet.specific_enthalpy_j_kg).unwrap();
            let geometry = input.droplet.geometry(input.law, t).unwrap();
            GasLiquidChamberSurface {
                id: input.id,
                law: input.law,
                binding: input.binding,
                liquid: fs_thermal::evaporation::LiquidSurfaceState {
                    mass_kg: input.droplet.mass_kg,
                    specific_enthalpy_j_kg: input.droplet.specific_enthalpy_j_kg,
                    area_m2: geometry.area_m2(),
                },
                temperature_from_enthalpy: input.temperature_from_enthalpy,
                mass_transfer_coefficient_m_s: geometry
                    .mass_transfer_coefficient_at_sherwood_m_s(
                        input.binary_diffusivity_m2_s,
                        input.sherwood_number,
                    )
                    .unwrap(),
            }
        })
        .collect();
    let direct =
        advance_chamber_with_gas(&original_surfaces, &original, &source, dt, solve(), || {
            false
        })
        .unwrap();
    assert_eq!(candidate.gas(), direct.gas());
    assert_eq!(candidate.step(), direct.step());
    let original_mass_defect = (direct.gas_exchange().step().mass_residual_kg.abs()
        + direct
            .surfaces()
            .iter()
            .map(|surface| surface.mass_residual_kg.abs())
            .sum::<f64>())
    .max(direct.step().mass_residual_kg.abs())
    .max(direct.step().maximum_species_mass_residual_kg.abs());
    near(
        candidate.report().sum_absolute_mass_residual_kg,
        original_mass_defect,
        2e-29,
    );
    assert!(candidate.report().sum_absolute_mass_residual_kg <= controls().gas.mass_tolerance_kg);
    let mut liquid_change = 0.0;
    let mut energy_change = 0.0;
    for (step, direct) in candidate.droplets().iter().zip(direct.surfaces()) {
        let input = inputs.iter().find(|input| input.id == step.id).unwrap();
        let t = temperature(input.droplet.specific_enthalpy_j_kg).unwrap();
        let p = original.state().state().pressure;
        let pv = original.species_partial_pressure_pa(&vapor()).unwrap();
        let ps = Phase.properties(t).unwrap().saturation_pressure_pa;
        // Independent integrated spherical Stefan rate: area*k = 2*pi*r*Sh*D.
        let request = dt
            * 2.0
            * core::f64::consts::PI
            * step.geometry[0].radius_m()
            * input.sherwood_number
            * input.binary_diffusivity_m2_s
            * (0.018015 * p / (MOLAR_GAS_CONSTANT * 320.0))
            * ((p - pv) / (p - ps)).ln();
        near(step.phase.requested_transferred_mass_kg, request, 2e-22);
        assert_eq!(&step.exchange, direct);
        assert_eq!(step.phase, step.exchange);
        assert_eq!(step.heat, None);
        assert_eq!(step.binary_diffusivity_m2_s, input.binary_diffusivity_m2_s);
        assert_eq!(step.sherwood_number, input.sherwood_number);
        assert_eq!(step.exchange.liquid.area_m2, step.geometry[0].area_m2());
        near(
            step.geometry[1].radius_m(),
            (3.0 * step.droplet.mass_kg / (4.0 * core::f64::consts::PI * 1000.0)).cbrt(),
            2e-17,
        );
        assert_eq!(
            step.geometry[1].radius_m() < step.geometry[0].radius_m(),
            step.phase.transferred_mass_kg > 0.0
        );
        let dm = step.droplet.mass_kg - input.droplet.mass_kg;
        liquid_change += dm;
        energy_change += step.droplet.mass_kg * step.droplet.specific_enthalpy_j_kg
            - input.droplet.mass_kg * input.droplet.specific_enthalpy_j_kg
            + binding.liquid_enthalpy_offset_j_kg * dm;
    }
    assert!(candidate.droplets()[0].phase.transferred_mass_kg > 0.0);
    assert!(candidate.droplets()[1].phase.transferred_mass_kg < 0.0);
    near(
        liquid_change + candidate.gas().mass_kg() - original.mass_kg(),
        0.0,
        controls().gas.mass_tolerance_kg,
    );
    near(
        liquid_change + candidate.gas().species_mass_kg(&vapor()).unwrap()
            - original.species_mass_kg(&vapor()).unwrap(),
        0.0,
        controls().gas.mass_tolerance_kg,
    );
    near(
        candidate.gas().species_mass_kg(&carrier()).unwrap(),
        original.species_mass_kg(&carrier()).unwrap(),
        controls().gas.mass_tolerance_kg,
    );
    energy_change += candidate.gas().mass_kg()
        * candidate.gas().state().specific_internal_energy_j_kg()
        - original.mass_kg() * original.state().specific_internal_energy_j_kg();
    near(energy_change, 0.0, controls().gas.energy_tolerance_j);
    assert!(candidate.step().sum_absolute_energy_residual_j <= controls().gas.energy_tolerance_j);
    assert_eq!(candidate.gas().volume_m3(), original.volume_m3());
    assert_eq!(candidate.gas().point().axes().get("local-lot"), Some(&11.0));
    assert!(std::ptr::eq(candidate.catalog(), &source));
    let next_inputs: Vec<_> = candidate
        .droplets()
        .iter()
        .map(|step| {
            let input = *inputs.iter().find(|input| input.id == step.id).unwrap();
            GasDropletSurface {
                droplet: step.droplet,
                ..input
            }
        })
        .collect();
    let next =
        advance_droplets_with_gas(&next_inputs, candidate.gas(), &source, dt, solve(), || {
            false
        })
        .unwrap();
    for (before, after) in candidate.droplets().iter().zip(next.droplets()) {
        assert_eq!(after.geometry[0], before.geometry[1]);
        assert_ne!(
            after.mass_transfer_coefficient_m_s,
            before.mass_transfer_coefficient_m_s
        );
    }
    assert_eq!(original, saved);
}

#[test]
fn singleton_retains_scalar_finite_kinetics_and_heat_and_empty_retains_gas() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let original = gas(&source, 330.0, 5000.0, 1e-4);
    let kinetic = KineticEvaporation::try_new(
        &Phase,
        InterfaceKinetics::HertzKnudsen {
            accommodation_coefficient: 1e-4,
        },
    )
    .unwrap();
    let input = GasDropletSurface {
        law: &kinetic,
        ..row(7, &binding, 1e-6, 300.0, D, 6.5)
    };
    for heated in [false, true] {
        let (candidate, scalar) = if heated {
            (
                advance_droplets_with_gas_and_heat(
                    &[GasDropletHeatSurface {
                        surface: input,
                        heat_transfer_coefficient_w_m2_k: 80.0,
                        heat_controls: heat_controls(),
                    }],
                    &original,
                    &source,
                    0.05,
                    solve(),
                    outer(),
                    || false,
                )
                .unwrap(),
                advance_droplet_with_gas_and_heat_at_sherwood(
                    &kinetic,
                    &binding,
                    input.droplet,
                    &original,
                    &source,
                    temperature,
                    D,
                    6.5,
                    80.0,
                    0.05,
                    controls(),
                    heat_controls(),
                    || false,
                )
                .unwrap(),
            )
        } else {
            (
                advance_droplets_with_gas(&[input], &original, &source, 0.05, solve(), || false)
                    .unwrap(),
                advance_droplet_with_gas_at_sherwood(
                    &kinetic,
                    &binding,
                    input.droplet,
                    &original,
                    &source,
                    temperature,
                    D,
                    6.5,
                    0.05,
                    controls(),
                    || false,
                )
                .unwrap(),
            )
        };
        let step = &candidate.droplets()[0];
        assert_eq!(step.droplet, scalar.droplet());
        assert_eq!(step.geometry, *scalar.geometry());
        assert_eq!(candidate.gas(), scalar.gas());
        assert_eq!(
            step.phase.transferred_mass_kg,
            scalar.phase_step().transferred_mass_kg
        );
        assert_eq!(
            step.phase.surface_temperature_k,
            scalar.phase_step().surface_temperature_k
        );
        assert_eq!(
            step.exchange.surface_temperature_k,
            scalar.step().surface_temperature_k
        );
        assert_eq!(
            candidate.step().energy_into_gas_j,
            scalar.step().energy_into_gas_j
        );
        assert_eq!(
            candidate.heat().map(|h| h.heat_into_liquid_j),
            scalar.heat().map(|h| h.heat_into_liquid_j)
        );
    }
    let empty =
        advance_droplets_with_gas(&[], &original, &source, 0.05, solve(), || false).unwrap();
    assert!(empty.droplets().is_empty());
    assert_eq!(empty.gas(), &original);
    assert_eq!(empty.step().transferred_mass_kg, 0.0);
    assert_eq!(empty.step().sum_absolute_energy_residual_j, 0.0);
}

struct ExpandingPhase;
impl EvaporationLaw for ExpandingPhase {
    fn temperature_range_k(&self) -> [f64; 2] {
        Phase.temperature_range_k()
    }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        let mut value = Phase.properties(t)?;
        value.liquid_density_kg_m3 = 1000.0 / (1.0 + 0.001 * (t - 300.0));
        Ok(value)
    }
}

#[test]
fn collective_heat_changes_each_actual_density_and_size_even_when_diffusion_is_sealed() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let original = gas(&source, 320.0, 1000.0, 1e-4);
    let inputs = [
        row(8, &binding, 1e-6, 280.0, 0.0, 2.0),
        row(3, &binding, 8e-6, 300.0, 0.0, 5.0),
        row(5, &binding, 27e-6, 340.0, 0.0, 8.0),
    ]
    .map(|surface| GasDropletHeatSurface {
        surface: GasDropletSurface {
            law: &ExpandingPhase,
            ..surface
        },
        heat_transfer_coefficient_w_m2_k: 80.0,
        heat_controls: heat_controls(),
    });
    let candidate = advance_droplets_with_gas_and_heat(
        &inputs,
        &original,
        &source,
        0.05,
        solve(),
        outer(),
        || false,
    )
    .unwrap();
    let reversed = advance_droplets_with_gas_and_heat(
        &[inputs[2], inputs[1], inputs[0]],
        &original,
        &source,
        0.05,
        solve(),
        outer(),
        || false,
    )
    .unwrap();
    assert_eq!(candidate, reversed);
    assert_eq!(candidate.gas().mass_kg(), original.mass_kg());
    assert_eq!(candidate.step().transferred_mass_kg, 0.0);
    assert_eq!(candidate.step().energy_into_gas_j, 0.0);
    let mut heat_sum = 0.0;
    for step in candidate.droplets() {
        let input = inputs
            .iter()
            .find(|input| input.surface.id == step.id)
            .unwrap()
            .surface;
        let heat = step.heat.unwrap();
        assert_eq!(step.droplet.mass_kg, input.droplet.mass_kg);
        assert_eq!(step.mass_transfer_coefficient_m_s, 0.0);
        assert_eq!(step.phase.surface_temperature_k[0], heat.temperature_k[1]);
        assert_eq!(step.phase.surface_temperature_k[1], heat.temperature_k[1]);
        assert_ne!(step.geometry[1].radius_m(), step.geometry[0].radius_m());
        let density = 1000.0 / (1.0 + 0.001 * (heat.temperature_k[1] - 300.0));
        assert_eq!(step.geometry[1].liquid_density_kg_m3(), density);
        near(
            step.geometry[1].radius_m(),
            (3.0 * step.droplet.mass_kg / (4.0 * core::f64::consts::PI * density)).cbrt(),
            2e-17,
        );
        near(
            step.exchange.liquid_enthalpy_change_j,
            input.droplet.mass_kg
                * (step.droplet.specific_enthalpy_j_kg - input.droplet.specific_enthalpy_j_kg),
            controls().gas.energy_tolerance_j,
        );
        heat_sum += heat.heat_into_liquid_j;
    }
    near(
        heat_sum,
        candidate.heat().unwrap().heat_into_liquid_j,
        controls().gas.energy_tolerance_j,
    );
    near(
        candidate.step().liquid_enthalpy_change_j + candidate.step().gas_internal_energy_change_j,
        0.0,
        controls().gas.energy_tolerance_j,
    );
    assert_eq!(candidate.phase_step().gas_internal_energy_change_j, 0.0);
}

#[test]
fn individually_admissible_condensation_shares_the_original_vapor_donor_once() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let original = gas(&source, 300.0, 5000.0, 1e-6);
    let saved = original.clone();
    let first = row(1, &binding, 8e-6, 300.0, D, 2.0);
    let second = GasDropletSurface { id: 2, ..first };
    let r = (3.0 * first.droplet.mass_kg / (4.0 * core::f64::consts::PI * 1000.0)).cbrt();
    let rate = 4.0
        * core::f64::consts::PI
        * r
        * D
        * (0.018015 * 100000.0 / (MOLAR_GAS_CONSTANT * 300.0))
        * ((100000.0_f64 - 5000.0) / (100000.0 - 1500.0)).ln();
    let dt = 0.2 * original.species_mass_kg(&vapor()).unwrap() / -rate;
    for input in [first, second] {
        assert!(
            advance_droplets_with_gas(&[input], &original, &source, dt, solve(), || false,).is_ok()
        );
    }
    let error =
        advance_droplets_with_gas(&[second, first], &original, &source, dt, solve(), || false)
            .unwrap_err();
    assert!(
        matches!(
            error,
            GasEvaporationError::StepSizeLimit {
                reason: "shared initial vapor inventory withdrawal limit",
                ..
            }
        ),
        "{error:?}"
    );
    assert_eq!(original, saved);
    assert!(
        advance_droplets_with_gas(
            &[first, second],
            &original,
            &source,
            dt / 4.0,
            solve(),
            || false,
        )
        .is_ok()
    );
}

#[test]
fn invalid_seeds_transport_ids_work_and_late_cancellation_do_not_publish_partial_spheres() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let original = gas(&source, 320.0, 5000.0, 1e-4);
    let saved = original.clone();
    let input = row(1, &binding, 1e-6, 300.0, D, 2.0);
    for bad in [
        GasDropletSurface {
            binary_diffusivity_m2_s: -D,
            ..input
        },
        GasDropletSurface {
            sherwood_number: 0.0,
            ..input
        },
        GasDropletSurface {
            droplet: SphericalDropletState {
                mass_kg: 0.0,
                ..input.droplet
            },
            ..input
        },
    ] {
        assert!(
            advance_droplets_with_gas(&[bad], &original, &source, 0.05, solve(), || false,)
                .is_err()
        );
    }
    assert!(
        advance_droplets_with_gas(&[input, input], &original, &source, 0.05, solve(), || false,)
            .is_err()
    );
    let second = GasDropletSurface { id: 2, ..input };
    assert!(
        advance_droplets_with_gas(
            &[input, second],
            &original,
            &source,
            0.05,
            GasLiquidChamberSolve {
                maximum_surfaces: 1,
                ..solve()
            },
            || false,
        )
        .is_err()
    );
    let mut calls = 0;
    advance_droplets_with_gas(&[input, second], &original, &source, 0.05, solve(), || {
        calls += 1;
        false
    })
    .unwrap();
    let last_gate = calls;
    let mut calls = 0;
    assert_eq!(
        advance_droplets_with_gas(&[input, second], &original, &source, 0.05, solve(), || {
            calls += 1;
            calls == last_gate
        },)
        .unwrap_err(),
        GasEvaporationError::Cancelled
    );
    assert_eq!(original, saved);
    assert!(
        advance_droplets_with_gas(&[input, second], &original, &source, 0.05, solve(), || {
            false
        },)
        .is_ok()
    );
}
