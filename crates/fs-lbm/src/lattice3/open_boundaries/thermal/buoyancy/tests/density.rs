use super::*;

#[derive(Clone, Copy)]
struct QuadraticDensity;
impl DensityDepartureLaw3 for QuadraticDensity {
    type Error = core::convert::Infallible;
    fn relative_density_departure(&self, theta: f64) -> Result<f64, Self::Error> {
        // Manufactured density maximum, not an empirical fluid fit.
        Ok(-0.04 * (theta - 0.5).powi(2))
    }
}

fn nonlinear_impulse<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let law =
        BoussinesqDensityLaw3::try_new(REFERENCE_RHO, GRAVITY, 0.01, QuadraticDensity).unwrap();
    let mut owner = make_with_policy::<L, Q, _>(
        collision,
        [4, 3, 2],
        BoundarySpec3::periodic(),
        [ThermalOpenFace3::Periodic; 6],
        0.25,
        law,
    );
    let cells = owner.thermal().populations().len();
    for cell in 0..cells {
        let (density, velocity) = owner.flow().macroscopic(cell).unwrap();
        near(density, RHO, 2e-14);
        for axis in 0..3 {
            near(velocity[axis], INITIAL_VELOCITY[axis], 2e-14);
        }
    }
    let heat = 0.0625;
    let mut expected_velocity = INITIAL_VELOCITY;
    for n in 0..8 {
        let opening = 0.25 + n as f64 * heat;
        let ending = opening + heat;
        // Independent original density relation and native impulse quadrature.
        // The thermal anomaly changes sign in slope as it crosses its maximum.
        let rho_delta_0 = -REFERENCE_RHO * 0.04 * (opening - 0.5).powi(2);
        let rho_delta_1 = -REFERENCE_RHO * 0.04 * (ending - 0.5).powi(2);
        let forces_0 =
            core::array::from_fn::<_, 3, _>(|axis| BACKGROUND[axis] + rho_delta_0 * GRAVITY[axis]);
        let forces_1 =
            core::array::from_fn::<_, 3, _>(|axis| BACKGROUND[axis] + rho_delta_1 * GRAVITY[axis]);
        let report = advance(&mut owner, &vec![heat; cells]);
        for axis in 0..3 {
            let impulse = 0.5 * forces_0[axis] + 0.5 * forces_1[axis];
            expected_velocity[axis] += impulse / RHO;
            near(report.body_impulse[axis], cells as f64 * impulse, 3e-14);
            near(
                report.background_impulse[axis],
                cells as f64 * BACKGROUND[axis],
                3e-15,
            );
        }
        near(report.heat.stored_change, cells as f64 * heat, 3e-13);
        near(report.heat.combined_defect, 0.0, 3e-13);
        for cell in 0..cells {
            near(owner.thermal().temperature(cell).unwrap(), ending, 3e-14);
            let (density, velocity) = owner.flow().macroscopic(cell).unwrap();
            near(density, RHO, 2e-13);
            for axis in 0..3 {
                near(velocity[axis], expected_velocity[axis], 3e-13);
                near(owner.flow().forces()[cell][axis], forces_1[axis], 3e-15);
            }
        }
        if n == 3 {
            assert!(report.buoyancy_impulse[1] > 0.0);
            near(owner.flow().forces()[0][1], BACKGROUND[1], 1e-25);
        }
    }
    near(owner.force_policy().reference_density(), REFERENCE_RHO, 0.0);
}

#[test]
fn nonlinear_density_crosses_its_maximum_with_original_initial_frame_and_trapezoidal_impulse() {
    nonlinear_impulse::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    nonlinear_impulse::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Trt {
            tau: 0.9,
            magic: 3.0 / 16.0,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    });
}

#[derive(Clone, Copy)]
struct InteriorDensityPeak;
impl DensityDepartureLaw3 for InteriorDensityPeak {
    type Error = core::convert::Infallible;
    fn relative_density_departure(&self, theta: f64) -> Result<f64, Self::Error> {
        Ok(0.02 / (1.0 + 40000.0 * (theta - 0.5).powi(2)))
    }
}

fn peak_fixture(
    maximum: f64,
) -> BoussinesqOpenGrid3<D3q19, 19, BoussinesqDensityLaw3<InteriorDensityPeak>> {
    let native = BoundarySpec3::new([
        FaceBoundary3::Periodic,
        FaceBoundary3::Periodic,
        FaceBoundary3::stationary_wall(),
        FaceBoundary3::stationary_wall(),
        FaceBoundary3::Periodic,
        FaceBoundary3::Periodic,
    ]);
    let faces = [
        ThermalOpenFace3::Periodic,
        ThermalOpenFace3::Periodic,
        ThermalOpenFace3::Wall(ThermalWall::HeatFlux(0.1)),
        ThermalOpenFace3::Wall(ThermalWall::Adiabatic),
        ThermalOpenFace3::Periodic,
        ThermalOpenFace3::Periodic,
    ];
    let initial: Vec<_> = (0..24).map(|cell| 0.1 * (cell % 4 % 3) as f64).collect();
    make_with_initial_temperatures(
        CollisionModel3::Bgk { tau: 0.9 },
        [4, 3, 2],
        native,
        faces,
        &initial,
        BoussinesqDensityLaw3::try_new(REFERENCE_RHO, GRAVITY, maximum, InteriorDensityPeak)
            .unwrap(),
    )
}

#[test]
fn nonlinear_admission_visits_interior_face_traces_and_actual_source_half_kicks() {
    let mut actual = peak_fixture(0.01);
    let mut witness = peak_fixture(0.03);
    let zero = vec![0.0; 24];
    let initial = accepted(&actual);
    let admitted = advance(&mut witness, &zero);
    let wall = admitted.heat.walls[2].unwrap();
    near(wall.min_temperature, 0.4, 3e-14);
    near(wall.max_temperature, 0.6, 3e-14);
    // An extrema-only law gate would admit BOTH endpoints of this same face.
    assert!(
        actual
            .force_policy()
            .sample(wall.min_temperature, [0.0; 3])
            .is_ok()
    );
    assert!(
        actual
            .force_policy()
            .sample(wall.max_temperature, [0.0; 3])
            .is_ok()
    );
    let refused = actual.step_with_poll(None, &zero, None, [0.0; 3], &mut || Ok(()));
    assert!(matches!(
        refused,
        Err(BoussinesqOpenError3::BoundaryLaw {
            face: 2,
            source: BoussinesqDensityLawError3::Law(BoussinesqLawError3::DensityDeparture { .. }),
        })
    ));
    assert_eq!(accepted(&actual), initial);

    // The source's first half kick reaches theta=0.5 even though its opening
    // and full-kick endpoint (theta=0 and 1) satisfy this original density law.
    let mut spike = zero.clone();
    spike[0] = 1.0;
    let refused = actual.step_with_poll(None, &spike, None, [0.0; 3], &mut || Ok(()));
    assert!(matches!(
        refused,
        Err(BoussinesqOpenError3::Law {
            cell: 0,
            source: BoussinesqDensityLawError3::Law(BoussinesqLawError3::DensityDeparture { .. }),
        })
    ));
    assert_eq!(accepted(&actual), initial);

    // Correcting the authored heat flux reuses original accepted scalar,
    // force and wall history; refused scratch cannot contaminate the retry.
    let mut reference = peak_fixture(0.01);
    for owner in [&mut actual, &mut reference] {
        owner
            .set_wall_boundary(2, ThermalWall::HeatFlux(0.05))
            .unwrap();
    }
    let expected = advance(&mut reference, &zero);
    let before = accepted(&actual);
    let mut polls = 0;
    let mut control = peak_fixture(0.01);
    control
        .set_wall_boundary(2, ThermalWall::HeatFlux(0.05))
        .unwrap();
    control
        .step_with_poll(None, &zero, None, [0.0; 3], &mut || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut calls = 0;
    let cancelled = actual.step_with_poll(None, &zero, None, [0.0; 3], &mut || {
        calls += 1;
        if calls == polls {
            Err(GridError3::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        cancelled,
        Err(BoussinesqOpenError3::Flow(GridError3::Cancelled))
    ));
    assert_eq!(accepted(&actual), before);
    assert_eq!(advance(&mut actual, &zero), expected);
    assert_eq!(accepted(&actual), accepted(&reference));
    near(expected.heat.stored_change, 8.0 * 0.05, 3e-14);
    near(expected.heat.combined_defect, 0.0, 3e-14);
}
