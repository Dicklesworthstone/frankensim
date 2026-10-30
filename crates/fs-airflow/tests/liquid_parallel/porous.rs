//! Original skeleton changes redistribute flow through the same network.
use super::*;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::porous::{
    LiquidResistanceLaw, PorousLiquidResistance,
};
use fs_material::porous::permeability::{
    DarcyGasBarrier, IntrinsicPermeability, PorousGasLayer, PorousGasSupport,
};

#[test]
fn porous_data_reaches_branch_split_heat_and_common_outlet_enthalpy() {
    with_source(false, false, |cx, law, source| {
        let a = path(cx, law, source, "porous-wall");
        let b = path(cx, law, source, "channel-wall");
        for factor in [1e10, 2e10] {
            let layer = PorousGasLayer::try_new(
                IntrinsicPermeability::isotropic(0.01 / (factor * geometry().flow_area.value()))
                    .unwrap(),
                0.01,
                0.0,
                "synthetic intrinsic pore geometry, not measurements",
            )
            .unwrap();
            let barrier = DarcyGasBarrier::try_new(
                vec![layer],
                [1.0, 0.0, 0.0],
                PorousGasSupport {
                    temperature_k: (250.0, 450.0),
                    pressure_pa: (90000.0, 1e6),
                    maximum_superficial_speed_m_s: 1.0,
                },
            )
            .unwrap();
            let resistance = PorousLiquidResistance::try_new(
                barrier.continuum_geometry().unwrap(),
                geometry().flow_area,
                Length::new(0.01),
                (0.0, 2300.0),
                Pressure::new(100000.0),
            )
            .unwrap();
            let ar = [LiquidResistanceSection {
                region: "porous-wall",
                resistance: LiquidResistanceLaw::Porous(resistance),
            }];
            let br = [LiquidResistanceSection {
                region: "channel-wall",
                resistance: LiquidResistanceLaw::Geometry(loss("channel-wall", 1.0)[0].resistance),
            }];
            let network = ParallelLiquidNetwork::try_new(
                cx,
                &[
                    ParallelLiquidBranch {
                        name: "porous",
                        path: &a,
                        resistances: &ar,
                        pressure_budget: branch_budget(),
                    },
                    ParallelLiquidBranch {
                        name: "channel",
                        path: &b,
                        resistances: &br,
                        pressure_budget: branch_budget(),
                    },
                ],
            )
            .unwrap();
            let result = network
                .solve_pressure_driven(cx, &[400.0, 280.0], Pressure::new(50.0), local(), budget())
                .unwrap();
            let masses = [50.0 * 1000.0 / (factor * 0.002), 0.0025];
            let mut total_heat = 0.0;
            for ((branch, mass), wall) in result.branches.iter().zip(masses).zip([400.0, 280.0]) {
                near(branch.solution.mass_flow_kg_per_s, mass);
                let outlet = wall + (300.0 - wall) * (-183.0 * 0.003 / (mass * 2000.0)).exp();
                near(
                    branch.solution.thermal.regions[0].outlet.temperature_k(),
                    outlet,
                );
                total_heat += mass * 2000.0 * (outlet - 300.0);
            }
            near(result.total_heat_rate_w, total_heat);
            near(
                result.mixed_outlet.temperature_k(),
                300.0 + total_heat / ((masses[0] + masses[1]) * 2000.0),
            );
            assert_eq!(
                result.branches[0].solution.hydraulics.regions[0]
                    .porous
                    .as_ref()
                    .unwrap()
                    .barrier,
                barrier
            );
            assert!(
                result.branches[1].solution.hydraulics.regions[0]
                    .porous
                    .is_none()
            );
        }
    });
}
