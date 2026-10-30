//! Manufactured melt-driven convection, checked against a native one-cycle balance.
use super::*;
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::phase::transport::{
    PhaseTransportSourceQuery, resolve_phase_transport_state_with_density_departure,
};
use fs_material::thermal::source::ThermalMaterialPoint;

fn density_source(ramp: f64, mu: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, 1.0),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            2.0,
        ),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, mu),
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
                abscissa: "liquid".into(),
                abscissa_dims: Dims([0; 6]),
                knots: vec![(0.0, 1.0 + ramp), (1.0, 1.0)],
                dims,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    card(claims)
}

fn buoyancy() -> BoussinesqConfig<3> {
    BoussinesqConfig {
        gravity_m_s2: [0.0, -10.0, 0.0],
        max_relative_density_departure: 0.025,
        max_relative_specific_heat_departure: 0.0,
    }
}

fn phase_source(
    card: &MaterialCard,
    law: &MushyPhaseLaw<'_>,
    h: f64,
) -> Result<fs_material::phase::transport::ResolvedPhaseTransportState, fs_material::MaterialError>
{
    let context = point(500.0, 100000.0).with("liquid", 0.0).unwrap();
    resolve_phase_transport_state_with_density_departure(
        law,
        h,
        PhaseTransportSourceQuery {
            source: ThermalMaterialPoint {
                material: card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: Some("liquid"),
            },
            viscosity_basis: LiquidViscosityBasis::Dynamic,
            extra_requirements: &[],
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
        0.025,
    )
}

#[derive(Clone)]
struct MeltRing;
impl DepositionObserver for MeltRing {
    type Error = DepositionError;
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        let target = [25.0, 50.0, 75.0, 25.0, -20.0, 75.0, 25.0, 50.0, 75.0];
        Ok(target
            .into_iter()
            .enumerate()
            .map(|(c, h)| {
                let watts = if thermal.time_s() < 0.1 {
                    (h + 20.0) / 0.1
                } else if (0.2..0.3).contains(&thermal.time_s()) && c == 4 {
                    2000.0
                } else {
                    0.0
                };
                (CellId(c as u64), watts)
            })
            .collect())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        Ok(())
    }
}

#[test]
fn cold_material_melts_then_drives_native_convection_and_opens_the_retained_center() {
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new_clock_free();
    let config = closed();
    let mut responses = Vec::new();
    for ramp in [0.0, 0.02] {
        let card = density_source(ramp, 1.0);
        let curve = chart(&card);
        let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
        let original = geometry(&domain, &curve, [-20.0; 9]);
        if ramp != 0.0 {
            assert!(
                MaterialFlowFeedback::try_new_phase(
                    &original,
                    &domain,
                    input(&ids, &curve),
                    &config,
                    (),
                    &gate,
                    10,
                    &law,
                    |_, cell: &DepositedCell, _| phase_source(&card, &law, cell.specific_h_j_kg),
                )
                .is_err()
            );
        }
        let history = MaterialFlowFeedback::try_new_phase_boussinesq(
            &original,
            &domain,
            input(&ids, &curve),
            &config,
            MeltRing,
            &gate,
            10,
            &law,
            FluidTransportState::try_new(1.0, 2.0, 0.1, 1.0).unwrap(),
            buoyancy(),
            |_, cell: &DepositedCell, _| phase_source(&card, &law, cell.specific_h_j_kg),
        )
        .unwrap()
        .with_unsteady_viscous_heating(controls())
        .unwrap();
        assert!(
            history
                .endpoint()
                .iter()
                .all(|s| s.material.transport_state().is_err())
        );
        let (mushy, history, _) = original
            .advance_observed(&history, 0.1, step(0.1), || Ok(()))
            .unwrap();
        assert_eq!(history.momentum_solves(), 0);
        assert_eq!(
            history.unsteady_heating_audit().unwrap().generated_heat_j,
            0.0
        );
        let cold = history.applied().unwrap();
        assert!(
            cold.buoyancy_force_n_m3
                .as_ref()
                .unwrap()
                .iter()
                .all(|f| *f == [0.0; 3])
        );
        for reaction in cold.buoyancy_support_reaction_n_m3.as_ref().unwrap() {
            assert!((reaction[1] - 10.0 * ramp).abs() < 1e-14);
        }
        let (moving, history, audit) = mushy
            .advance_observed(&history, 0.2, step(0.1), || Ok(()))
            .unwrap();
        let applied = history.applied().unwrap();
        let native = history.unsteady_flow().unwrap();
        // Eight MAC faces form one cycle. Full stress gives 48*mu*q²;
        // Darcy adds sum(beta_cell)*q², beta=mu*(1-phi)^2/phi^3.
        // Reference face inertia adds 8*q/dt; density contrast gives -10*ramp.
        let beta = |phi: f64| (1.0 - phi).powi(2) / phi.powi(3);
        let resistance = 48.0 + 3.0 * beta(0.25) + 2.0 * beta(0.5) + 3.0 * beta(0.75);
        let q = -10.0 * ramp / (8.0 / 0.1 + resistance);
        let speed = native.field.flux_y(0, 1, 0);
        assert!((speed - q).abs() < 1e-8, "native {speed}, independent {q}");
        responses.push(speed);
        assert_eq!(applied.buoyancy_force_n_m3.as_ref().unwrap()[4], [0.0; 3]);
        assert_eq!(native.velocity_m_s[4], [0.0; 3]);
        let mechanical = history.unsteady_heating_audit().unwrap();
        assert!((mechanical.supplied_buoyancy_work_j - 0.1 * (-10.0 * ramp) * q).abs() < 1e-9);
        assert!((mechanical.final_kinetic_energy_j - 4.0 * q * q).abs() < 1e-9);
        assert!((mechanical.backward_euler_loss_j - 4.0 * q * q).abs() < 1e-9);
        assert!((mechanical.generated_heat_j - 0.1 * resistance * q * q).abs() < 1e-9);
        assert!((moving.energy_j() - mushy.energy_j() - mechanical.generated_heat_j).abs() < 1e-8);
        assert!((audit.external_heat_j - mechanical.generated_heat_j).abs() < 1e-8);
        assert_eq!(mechanical.support_capture_heat_j, 0.0);
        assert_eq!(moving.mass_kg(), original.mass_kg());
        let (melted, history, _) = moving
            .advance_observed(&history, 0.3, step(0.1), || Ok(()))
            .unwrap();
        assert!(history.endpoint()[4].transport.mobile_transport().is_some());
        let (_, history, _) = melted
            .advance_observed(&history, 0.4, step(0.1), || Ok(()))
            .unwrap();
        assert_eq!(history.unsteady_time_s(), Some(0.4));
        if ramp != 0.0 {
            assert!(
                history.unsteady_flow().unwrap().velocity_m_s[4]
                    .iter()
                    .any(|u| u.abs() > 1e-10)
            );
        }
    }
    assert!(responses[0].abs() < 1e-12);
    assert!(responses[1] < -1e-5);
}

#[test]
fn endpoint_density_departure_refuses_without_advancing_heat_or_native_momentum() {
    let card = density_source(0.02, 1.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| Some(CellId(c as u64)));
    let config = closed();
    let gate = CancelGate::new_clock_free();
    let original = geometry(&domain, &curve, [150.0; 9]);
    let mut bounded = buoyancy();
    bounded.max_relative_density_departure = 0.001;
    let history = MaterialFlowFeedback::try_new_phase_boussinesq(
        &original,
        &domain,
        input(&ids, &curve),
        &config,
        CellPower {
            cell: CellId(4),
            watts: -1000.0,
        },
        &gate,
        10,
        &law,
        FluidTransportState::try_new(1.0, 2.0, 0.1, 1.0).unwrap(),
        bounded,
        |_, cell: &DepositedCell, _| phase_source(&card, &law, cell.specific_h_j_kg),
    )
    .unwrap()
    .with_unsteady_viscous_heating(controls())
    .unwrap();
    let before = (
        original.cells().clone(),
        history.unsteady_flow().unwrap().clone(),
    );
    assert!(
        original
            .advance_observed(&history, 0.1, step(0.1), || Ok(()))
            .is_err()
    );
    assert_eq!(original.cells(), &before.0);
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(history.unsteady_flow().unwrap(), &before.1);
    assert_eq!(history.unsteady_time_s(), Some(0.0));
    assert!(history.applied().is_none());
}
