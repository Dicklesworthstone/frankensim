//! Original finite reservoir/UV coupling. Trials allocate no product receipts.

use super::*;
use fs_math::dd::Dd;
use fs_qty::{Energy, Volume};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::rocket::equilibrium::blowdown::prepared::EquilibriumReservoirDepletion;
use fs_thermochem::equilibrium::composition::{
    isentrope::{EquilibriumIsentropeSampleV1, EquilibriumIsentropeV1},
    volume::energy_path::{EquilibriumEnergyVolumePathV1,
        anchored::EquilibriumAnchoredEnergyVolumeSampleV1},
};

pub(super) struct Trial {
    pub(super) donor: EquilibriumIsentropeSampleV1,
    pub(super) receiver: EquilibriumAnchoredEnergyVolumeSampleV1,
    outflow: EquilibriumBlowdownStep,
    remaining_mass_kg: f64,
    receiver_mass_kg: f64,
    receiver_increment: Dd,
    mass_residual: Dd,
    donor_enthalpy_change: f64,
}

pub(super) struct Workspace<'a> {
    pair: &'a EquilibriumGasPair,
    donor: usize,
    pub(super) depletion: EquilibriumReservoirDepletion<'a>,
    receiver: EquilibriumEnergyVolumePathV1<'a>,
    basis_mass: f64,
    donor_minus_receiver_h: f64,
    donor_minus_receiver_s: f64,
    receiver_flow_work: Dd,
    receiver_u: Dd,
    density_queries: usize,
    receiver_pressure_queries: usize,
    receiver_composition_iterations: usize,
    pub(super) nozzle_queries: usize,
}

fn finite(value: Dd, message: &'static str) -> Result<f64, Error> {
    let rounded = value.to_f64();
    if !value.hi.is_finite() || !value.lo.is_finite() || !rounded.is_finite()
        || ((value.hi != 0.0 || value.lo != 0.0) && rounded == 0.0) {
        Err(Error::Invalid(message))
    } else { Ok(rounded) }
}
fn flow_work(tank: &EquilibriumGasReservoir) -> Dd {
    Dd::from_f64(tank.isentropic_sample().pressure.value())
        * Dd::from_f64(tank.volume_m3()) / Dd::from_f64(tank.mass_kg())
}

impl<'a> Workspace<'a> {
    pub(super) fn new(
        pair: &'a EquilibriumGasPair, donor: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let d = &pair.tanks[donor];
        let r = &pair.tanks[1 - donor];
        let mut depletion = d.prepare_depletion(pair.controls.reservoir, &mut *cancelled)?;
        let path = depletion.isentrope();
        path.sample(d.isentropic_sample().pressure, &mut *cancelled)?;
        let basis_mass = path.feed_mass_kg();
        // Remove this re-query's own residual before comparing the two accepted
        // tank anchors. No absolute H or S references are subtracted.
        let dh = path.specific_enthalpy_change_from(r.caloric_anchor(), &mut *cancelled)?;
        let self_h = path.specific_enthalpy_change_from(d.caloric_anchor(), &mut *cancelled)?;
        let ds = path.specific_entropy_change_from(r.caloric_anchor(), &mut *cancelled)?;
        let self_s = path.specific_entropy_change_from(d.caloric_anchor(), &mut *cancelled)?;
        if (self_h * d.mass_kg()).abs() > pair.controls.reservoir.energy_tolerance_j
            || (self_s * d.mass_kg()).abs() > pair.controls.entropy_tolerance_j_k {
            return Err(Error::Invalid("re-prepared donor changed its accepted source state"));
        }
        let donor_minus_receiver_h = finite(Dd::from_f64(dh) - Dd::from_f64(self_h),
            "unrepresentable direct inter-tank enthalpy difference")?;
        let donor_minus_receiver_s = finite(Dd::from_f64(ds) - Dd::from_f64(self_s),
            "unrepresentable direct inter-tank entropy difference")?;
        let volume = finite(Dd::from_f64(r.volume_m3()) * Dd::from_f64(basis_mass)
            / Dd::from_f64(r.mass_kg()), "unrepresentable receiver feed-basis volume")?;
        let receiver = EquilibriumEnergyVolumePathV1::new(
            &pair.chemistry, r.isentropic_sample().temperature, Volume::new(volume),
            &pair.feed, pair.controls.receiver, &mut *cancelled,
        )?;
        if receiver.feed_mass_kg() != basis_mass {
            return Err(Error::Invalid("two-volume source normalization changed"));
        }
        let receiver_flow_work = flow_work(r);
        let receiver_u = Dd::from_f64(r.specific_enthalpy_j_kg()) - receiver_flow_work;
        let receiver_pressure_queries = receiver.reference_state().pressure_iterations();
        let receiver_composition_iterations = receiver.reference_state().equilibrium().composition_iterations();
        poll(cancelled)?;
        Ok(Self {
            pair, donor, depletion, receiver, basis_mass,
            donor_minus_receiver_h, donor_minus_receiver_s,
            receiver_flow_work, receiver_u, density_queries: 0,
            receiver_pressure_queries, receiver_composition_iterations, nozzle_queries: 0,
        })
    }

    pub(super) fn trial(
        &mut self, remaining_mass_kg: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial, Error> {
        poll(cancelled)?;
        let old = &self.pair.tanks[1 - self.donor];
        let outflow = self.depletion.sample_withdrawal(remaining_mass_kg, &mut *cancelled)?;
        self.density_queries += outflow.pressure_queries;
        let donor = self.depletion.isentrope().current_sample()
            .ok_or(Error::Invalid("finite donor lost its accepted sample"))?;
        let donor_enthalpy_change = self.depletion.isentrope().specific_enthalpy_change_from(
            self.pair.tanks[self.donor].caloric_anchor(), &mut *cancelled,
        )?;
        let dm = outflow.expelled_mass_kg;
        let receiver_mass_kg = old.mass_kg() + dm;
        if !positive(receiver_mass_kg) || receiver_mass_kg == old.mass_kg() {
            return Err(Error::Invalid("received mass is not representable"));
        }
        let mass_residual = Dd::from_f64(receiver_mass_kg)
            - Dd::from_f64(old.mass_kg()) - Dd::from_f64(dm);
        let mass_residual_kg = finite(mass_residual, "unrepresentable pair mass residual")?;
        if mass_residual_kg.abs() > self.pair.controls.reservoir.mass_tolerance_kg {
            return Err(Error::MassBalance { residual_kg: mass_residual_kg });
        }
        // U_r,new-U_r,old = dm*h_d + C, where C is the ORIGINAL finite
        // reservoir correction. Cancel the common formation datum first:
        // m_r,new*(u_r,new-u_r,old) = dm*(h_d-h_r+p_r/rho_r)+C.
        // Keep actual represented mass rounding and its datum sensitivity.
        let mass_reference = mass_residual * self.receiver_u;
        let reference_error = finite(mass_reference, "unrepresentable mass/reference energy")?;
        if reference_error.abs() > self.pair.controls.reservoir.energy_tolerance_j {
            return Err(Error::EnergyBalance { residual_j: reference_error });
        }
        let receiver_increment = Dd::from_f64(dm)
            * (Dd::from_f64(self.donor_minus_receiver_h) + self.receiver_flow_work)
            + Dd::from_f64(outflow.enthalpy_correction_j) - mass_reference;
        let target = finite(receiver_increment * Dd::from_f64(self.basis_mass)
            / Dd::from_f64(receiver_mass_kg), "unrepresentable receiver relative-energy target")?;
        let volume = finite(Dd::from_f64(old.volume_m3()) * Dd::from_f64(self.basis_mass)
            / Dd::from_f64(receiver_mass_kg), "unrepresentable received density")?;
        let receiver = self.receiver.sample_from_caloric_anchor(
            old.caloric_anchor(), Volume::new(volume), Energy::new(target), &mut *cancelled,
        )?;
        self.receiver_pressure_queries += receiver.sample.pressure_iterations;
        self.receiver_composition_iterations += receiver.sample.composition_iterations;
        let actual_residual = finite(
            Dd::from_f64(receiver.internal_energy_change_j)
                * Dd::from_f64(receiver_mass_kg) / Dd::from_f64(self.basis_mass)
                - receiver_increment,
            "unrepresentable physical receiving energy residual",
        )?;
        if actual_residual.abs() > self.pair.controls.reservoir.energy_tolerance_j {
            return Err(Error::EnergyBalance { residual_j: actual_residual });
        }
        poll(cancelled)?;
        Ok(Trial { donor, receiver, outflow, remaining_mass_kg,
            receiver_mass_kg, receiver_increment, mass_residual, donor_enthalpy_change })
    }

    // Exact thermodynamic tangent of the same quasistatic trajectory, used
    // only as a safeguarded proposal. Receiving entropy obeys
    // m*T*ds/dm = h_d-h_r; hence dp/dm = c^2/V + (alpha/kappa)/cv*(h_d-h_r)/m.
    // The ordinary subsonic throat responds to BOTH stagnation H and back p.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn rate_derivative(
        &self, trial: &Trial, mass_rate: f64, velocity: f64,
        exit_density: f64, exit_sound: f64, choked: bool,
    ) -> Option<f64> {
        let speed = if choked { exit_sound } else { velocity };
        let ratio = trial.donor.response.sound_speed_m_s / speed;
        let mut derivative = mass_rate / trial.remaining_mass_kg * ratio * ratio;
        if !choked {
            let r = &self.pair.tanks[1 - self.donor];
            let response = trial.receiver.sample.response;
            let flow_work = trial.receiver.sample.pressure.value() * r.volume_m3()
                / trial.receiver_mass_kg;
            let receiver_dh = trial.receiver.internal_energy_change_j / self.basis_mass
                + flow_work - self.receiver_flow_work.to_f64();
            let difference = self.donor_minus_receiver_h + trial.donor_enthalpy_change - receiver_dh;
            let cv = response.heat_capacity_v_j_k / self.basis_mass;
            let receiver_slope = response.sound_speed_m_s.powi(2) / r.volume_m3()
                + (response.thermal_expansion_per_k / response.isothermal_compressibility_per_pa)
                    / cv * difference / trial.receiver_mass_kg;
            derivative += mass_rate / (exit_density * velocity * velocity)
                * (1.0 - (velocity / exit_sound).powi(2)) * receiver_slope;
        }
        (derivative.is_finite() && derivative >= 0.0).then_some(derivative)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn finish(
        mut self, trial: Trial, midpoint_nozzle: EquilibriumNozzleFlow,
        duration_s: f64, quadrature_residual: f64, iterations: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EquilibriumGasPairCandidate, Error> {
        let pair = self.pair;
        let old_donor = &pair.tanks[self.donor];
        let old_receiver = &pair.tanks[1 - self.donor];
        let next = self.depletion.finish_remaining_mass(trial.remaining_mass_kg, &mut *cancelled)?;
        self.density_queries += next.step().pressure_queries;
        let (donor, outflow) = next.into_parts();
        // Reconstruct the actual receiving physical owner. Its new entropy
        // becomes its new birth reference; an old outflow isentrope is not reused.
        let receiver = EquilibriumGasReservoir::new(
            pair.chemistry.clone(), pair.feed.clone(), old_receiver.volume_m3(),
            trial.receiver.sample.temperature, trial.receiver.sample.pressure,
            pair.controls.nozzle.isentropic, &mut *cancelled,
        )?.with_conserved_mass(
            trial.receiver_mass_kg, pair.controls.reservoir.mass_tolerance_kg,
            &mut *cancelled,
        )?;
        let donor_pressure_pa = donor.isentropic_sample().pressure.value();
        let receiver_pressure_pa = receiver.isentropic_sample().pressure.value();
        if donor_pressure_pa < receiver_pressure_pa {
            return Err(Error::PressureCrossing { donor_pressure_pa, receiver_pressure_pa });
        }
        // The native UV tolerance cannot silently replace the actual returned
        // SP reservoir's energy/entropy/EOS admission.
        let receiver_path = EquilibriumIsentropeV1::new(
            &pair.chemistry, receiver.isentropic_sample().temperature,
            receiver.isentropic_sample().pressure, &pair.feed,
            pair.controls.nozzle.isentropic, &mut *cancelled,
        )?;
        self.receiver_composition_iterations += receiver.thermodynamics().equilibrium().composition_iterations()
            + receiver_path.reference_state().equilibrium().composition_iterations();
        if receiver_path.reference_state() != receiver.thermodynamics() {
            return Err(Error::Invalid("receiving publication changed its native TP source"));
        }
        let receiver_dh = receiver_path.specific_enthalpy_change_from(
            old_receiver.caloric_anchor(), &mut *cancelled,
        )?;
        let receiver_ds = receiver_path.specific_entropy_change_from(
            old_receiver.caloric_anchor(), &mut *cancelled,
        )?;
        let donor_ds = self.depletion.isentrope().specific_entropy_change_from(
            old_donor.caloric_anchor(), &mut *cancelled,
        )?;
        let receiver_du = Dd::from_f64(receiver_dh)
            - (flow_work(&receiver) - self.receiver_flow_work);
        // A final density reconstruction may have a different numerical warm
        // start. Its ORIGINAL finite correction, not the trial's, is decisive.
        let requested = trial.receiver_increment
            + (Dd::from_f64(outflow.enthalpy_correction_j)
                - Dd::from_f64(trial.outflow.enthalpy_correction_j));
        let energy_residual_j = finite(
            Dd::from_f64(receiver.mass_kg()) * receiver_du - requested,
            "unrepresentable actual two-volume energy residual",
        )?;
        if energy_residual_j.abs() > pair.controls.reservoir.energy_tolerance_j {
            return Err(Error::EnergyBalance { residual_j: energy_residual_j });
        }
        let entropy_production_j_k = finite(
            Dd::from_f64(donor.mass_kg()) * Dd::from_f64(donor_ds)
                + Dd::from_f64(receiver.mass_kg()) * Dd::from_f64(receiver_ds)
                - Dd::from_f64(outflow.expelled_mass_kg)
                    * Dd::from_f64(self.donor_minus_receiver_s)
                + trial.mass_residual * Dd::from_f64(
                    old_receiver.thermodynamics().entropy_j_k() / self.basis_mass),
            "unrepresentable actual pair entropy production",
        )?;
        if entropy_production_j_k < -pair.controls.entropy_tolerance_j_k {
            return Err(Error::Entropy { production_j_k: entropy_production_j_k });
        }
        let mut enthalpy_changes = [0.0; 2];
        let mut entropy_changes = [0.0; 2];
        enthalpy_changes[self.donor] = self.depletion.isentrope()
            .specific_enthalpy_change_from(pair.reference.tanks[self.donor].caloric_anchor(), &mut *cancelled)?;
        entropy_changes[self.donor] = self.depletion.isentrope()
            .specific_entropy_change_from(pair.reference.tanks[self.donor].caloric_anchor(), &mut *cancelled)?;
        enthalpy_changes[1 - self.donor] = receiver_path
            .specific_enthalpy_change_from(pair.reference.tanks[1 - self.donor].caloric_anchor(), &mut *cancelled)?;
        entropy_changes[1 - self.donor] = receiver_path
            .specific_entropy_change_from(pair.reference.tanks[1 - self.donor].caloric_anchor(), &mut *cancelled)?;
        let flow = resolve_on_isentrope(
            self.depletion.isentrope(), donor.isentropic_sample(), donor.caloric_anchor(),
            receiver_pressure_pa, midpoint_nozzle.effective_area.value(),
            pair.controls.nozzle, &mut *cancelled,
        )?;
        self.nozzle_queries += flow.work.pressure;
        let mut tanks = pair.tanks.clone();
        tanks[self.donor] = donor;
        tanks[1 - self.donor] = receiver;
        let element_residual = elements(pair, &tanks, cancelled)?;
        let reference_balance = reference_balances(pair, &tanks, enthalpy_changes, entropy_changes, self.basis_mass)?;
        let mass_residual_kg = finite(trial.mass_residual, "unrepresentable final pair mass residual")?;
        poll(cancelled)?;
        Ok(EquilibriumGasPairCandidate {
            endpoint: EquilibriumGasPair {
                chemistry: pair.chemistry.clone(), feed: pair.feed.clone(),
                tanks, reference: pair.reference.clone(), reference_balance,
                reference_element_residual: element_residual, controls: pair.controls,
            },
            step: EquilibriumGasPairStep {
                duration_s,
                direction: Some(if self.donor == 0 { EquilibriumGasPairDirection::LeftToRight }
                    else { EquilibriumGasPairDirection::RightToLeft }),
                donor: Some(outflow), midpoint_nozzle: Some(midpoint_nozzle),
                mass_residual_kg, maximum_relative_element_residual: element_residual,
                energy_residual_j, entropy_production_j_k,
                mass_change_since_reference_kg: reference_balance[0],
                energy_change_since_reference_j: reference_balance[1],
                entropy_change_since_reference_j_k: reference_balance[2],
                mass_quadrature_residual_kg: quadrature_residual, mass_iterations: iterations,
                density_pressure_queries: self.density_queries,
                receiver_pressure_queries: self.receiver_pressure_queries,
                receiver_composition_iterations: self.receiver_composition_iterations,
                nozzle_pressure_queries: self.nozzle_queries,
            },
        })
    }
}

fn elements(
    pair: &EquilibriumGasPair, tanks: &[EquilibriumGasReservoir; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let matrix = pair.chemistry.elemental_matrix();
    let mut maximum = 0.0_f64;
    for row in 0..matrix.row_count() {
        poll(cancelled)?;
        let count = |set: &[EquilibriumGasReservoir; 2]| {
            let mut value = Dd::ZERO;
            for tank in set {
                for (i, species) in tank.thermodynamics().equilibrium().species().iter().enumerate() {
                    value = value + Dd::from_f64(species.amount.value())
                        * Dd::from_f64(matrix.get(row, i).expect("original element matrix") as f64)
                        * Dd::from_f64(tank.physical_amount_scale());
                }
            }
            value
        };
        let before = count(&pair.reference.tanks);
        let after = count(tanks);
        let residual = finite(after - before, "unrepresentable actual elemental transfer")?;
        let relative = if before.to_f64() == 0.0 {
            if residual != 0.0 { return Err(Error::Invalid("absent elemental inventory was created")); }
            0.0
        } else {
            finite((after - before) / before, "unrepresentable relative elemental transfer")?.abs()
        };
        maximum = maximum.max(relative);
    }
    if maximum > pair.controls.relative_element_tolerance {
        return Err(Error::ElementBalance { maximum_relative_residual: maximum });
    }
    Ok(maximum)
}

fn reference_balances(
    pair: &EquilibriumGasPair, tanks: &[EquilibriumGasReservoir; 2],
    dh: [f64; 2], ds: [f64; 2], basis_mass: f64,
) -> Result<[f64; 3], Error> {
    let birth = &pair.reference.tanks;
    let dm = [0, 1].map(|i| Dd::from_f64(tanks[i].mass_kg()) - Dd::from_f64(birth[i].mass_kg()));
    let mass_change = dm[0] + dm[1];
    let mut energy = dm[1] * Dd::from_f64(pair.reference.right_minus_left_h)
        + mass_change * Dd::from_f64(birth[0].specific_enthalpy_j_kg());
    let mut entropy = dm[1] * Dd::from_f64(pair.reference.right_minus_left_s)
        + mass_change * Dd::from_f64(birth[0].thermodynamics().entropy_j_k() / basis_mass);
    for i in 0..2 {
        energy = energy + Dd::from_f64(tanks[i].mass_kg()) * Dd::from_f64(dh[i])
            - Dd::from_f64(tanks[i].volume_m3())
                * (Dd::from_f64(tanks[i].isentropic_sample().pressure.value())
                    - Dd::from_f64(birth[i].isentropic_sample().pressure.value()));
        entropy = entropy + Dd::from_f64(tanks[i].mass_kg()) * Dd::from_f64(ds[i]);
    }
    let mass_kg = finite(mass_change, "unrepresentable pair reference mass balance")?;
    let energy_j = finite(energy, "unrepresentable pair reference energy balance")?;
    let entropy_j_k = finite(entropy, "unrepresentable pair reference entropy balance")?;
    if mass_kg.abs() > pair.controls.reservoir.mass_tolerance_kg {
        return Err(Error::MassBalance { residual_kg: mass_kg });
    }
    if energy_j.abs() > pair.controls.reservoir.energy_tolerance_j {
        return Err(Error::EnergyBalance { residual_j: energy_j });
    }
    if entropy_j_k < -pair.controls.entropy_tolerance_j_k {
        return Err(Error::Entropy { production_j_k: entropy_j_k });
    }
    Ok([mass_kg, energy_j, entropy_j_k])
}
