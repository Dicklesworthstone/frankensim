//! Data-only material substitution reaches the existing reacting pipe network.
use super::*;
use crate::multicomponent::network::MixtureNetwork1d;
use crate::multicomponent::network::reactive::{ReactingNetwork1d, ReactiveNetworkError};
use crate::multicomponent::network::wall_exchange::{CircularLaminarWall, NetworkWallOperator};
use crate::multicomponent::wall_exchange::{
    WallContact, WallControl, WallDrag, WallExchangeError, WallThermal,
};

fn geometry() -> CircularLaminarWall {
    CircularLaminarWall {
        diameter: Length::new((0.7 / core::f64::consts::FRAC_PI_4).sqrt()),
        nusselt_number: 3.66,
        reynolds_number: [0.0, 2300.0],
        maximum_mach_number: 0.3,
    }
}

fn wall(n: &MixtureNetwork1d, material: ThermochemicalGasModel) -> NetworkWallOperator {
    NetworkWallOperator::new(
        n,
        vec![vec![
            WallContact {
                // These old authored coefficients are intentionally unrelated to
                // either source. The bound path must replace them at each half.
                drag: WallDrag::Linear { rate_per_s: 999.0 },
                thermal: WallThermal::Lumped {
                    temperature_k: 450.0,
                    heat_capacity_j_per_k: 50.0,
                    conductance_w_per_k: 999.0,
                },
            };
            n.branches()[0].cell_count()
        ]],
        Vec::new(),
        || false,
    )
    .unwrap()
    .with_circular_laminar_transport(adapter(material), vec![Some(geometry())], || false)
    .unwrap()
}

#[test]
fn material_cards_and_reaction_products_refresh_wall_drag_and_heat_in_same_step() {
    let material = model(0.02, 5.0, 2000.0);
    let initial = flow(&material, 2, |_| (0.1, 400.0));
    let n = MixtureNetwork1d::new(vec![initial], Vec::new(), Vec::new(), || false).unwrap();
    let mut wall = wall(&n, material.clone());
    let (_, mechanism, _) = fixture(0.0, 2000.0);
    let before = n.totals();
    let initial_energy = before.absolute_energy(n.gas()).unwrap() + 2.0 * 50.0 * 450.0;
    let mut reactive =
        ReactingNetwork1d::new(n, mechanism, controls().chemistry, || false).unwrap();
    let report = reactive
        .step_with_wall(
            &mut wall,
            Time::new(1e-4),
            ReactiveStepControl::default(),
            WallControl::default(),
            || false,
        )
        .unwrap();
    let current = reactive.branches()[0].cell(0).unwrap();
    assert!(
        current.partial_densities[1] > 0.0,
        "the real reactor must create product"
    );
    let primitive = reactive.gas().primitive(current, || false).unwrap();
    // Both fixture species have the same molar mass, so their mass fractions
    // are exactly their mole fractions. This independent query checks the
    // coefficients used after chemistry against the shared original material.
    let x = [
        current.partial_densities[0] / primitive.density,
        current.partial_densities[1] / primitive.density,
    ];
    let mut reference = adapter(material);
    let mut binary = [0.0; 4];
    let coefficient = reference
        .evaluate(
            Temperature::new(primitive.temperature),
            Pressure::new(primitive.pressure),
            &x,
            &mut binary,
            &mut || false,
        )
        .unwrap();
    let diameter = geometry().diameter.value();
    let expected_rate =
        32.0 * coefficient.dynamic_viscosity.value() / primitive.density / diameter / diameter;
    let expected_conductance = 3.66
        * coefficient.thermal_conductivity.value()
        * (4.0 * reactive.branches()[0].mesh().volumes()[0] / diameter)
        / diameter;
    assert!(coefficient.dynamic_viscosity.value() > 0.02);
    assert!(coefficient.thermal_conductivity.value() > 5.0);
    let mut wall_energy = 0.0;
    for contact in wall.branch_contacts(0).unwrap() {
        match contact.drag {
            WallDrag::Linear { rate_per_s } => close(rate_per_s, expected_rate, 2e-10),
            _ => panic!("material pipe drag missing"),
        }
        match contact.thermal {
            WallThermal::Lumped {
                temperature_k,
                heat_capacity_j_per_k,
                conductance_w_per_k,
            } => {
                close(conductance_w_per_k, expected_conductance, 2e-10);
                wall_energy += temperature_k * heat_capacity_j_per_k;
            }
            _ => panic!("finite wall storage missing"),
        }
    }
    close(reactive.totals().mass(), before.mass(), 1e-12);
    close(
        reactive.totals().absolute_energy(reactive.gas()).unwrap() + wall_energy,
        initial_energy,
        2e-10,
    );
    close(report.absolute_energy_residual_j, 0.0, 1e-7);
    assert!(report.wall.is_some());
    assert_eq!(wall.last_report().unwrap().elapsed, report.elapsed);
}

#[test]
fn new_product_material_refusal_rolls_back_whole_reacting_network_and_finite_walls() {
    // The absent product is outside its material domain at 400 K. Exact zero
    // initially permits A; the positive product created by chemistry must then
    // query B and refuse, rather than retaining pre-reaction wall coefficients.
    let material = model(0.02, 5.0, 399.0);
    let initial = flow(&material, 2, |_| (0.1, 400.0));
    let before = cells(&initial);
    let n = MixtureNetwork1d::new(vec![initial], Vec::new(), Vec::new(), || false).unwrap();
    let mut wall = wall(&n, material);
    let contacts = wall.branch_contacts(0).unwrap().to_vec();
    let totals = n.totals();
    let (_, mechanism, _) = fixture(0.0, 2000.0);
    let mut reactive =
        ReactingNetwork1d::new(n, mechanism, controls().chemistry, || false).unwrap();
    let error = reactive
        .step_with_wall(
            &mut wall,
            Time::new(1e-4),
            ReactiveStepControl {
                maximum_coupling_retries: 0,
                ..ReactiveStepControl::default()
            },
            WallControl::default(),
            || false,
        )
        .unwrap_err();
    fn material_refusal(error: &ReactiveNetworkError) -> bool {
        match error {
            ReactiveNetworkError::Wall(WallExchangeError::Transport(
                DissipationError::Material(_),
            )) => true,
            ReactiveNetworkError::RetryLimit { last, .. } => material_refusal(last),
            _ => false,
        }
    }
    assert!(material_refusal(&error), "unexpected refusal: {error:?}");
    assert_eq!(cells(&reactive.branches()[0]), before);
    assert_eq!(reactive.totals(), totals);
    assert_eq!(reactive.time(), Time::new(0.0));
    assert_eq!(reactive.last_report(), None);
    assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
    assert_eq!(wall.last_report(), None);
}
