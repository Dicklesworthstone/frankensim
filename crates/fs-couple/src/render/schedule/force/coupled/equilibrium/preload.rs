//! Nonlinear contact-loaded equilibrium, using the existing normal-reaction solve.
use super::*;
use super::super::contact::{ModalContact, ModalContactConfig, contact_column, force_tolerance};
use super::super::contact::multiple::MultiContactConfig;
use super::normal::{admit_setup, solve_joint};
use fs_dcontact::ContactStorage;
use fs_phs::Storage;

struct ZeroStorage;
impl Storage for ZeroStorage {
    fn hamiltonian(&self, _: &[f64]) -> f64 { 0.0 }
    fn gradient(&self, _: &[f64], out: &mut [f64]) { out.fill(0.0); }
}
struct Point {
    column: Vec<f64>,
    storage: ContactStorage,
    config: ModalContactConfig,
    gap: f64,
}
impl Point {
    fn force(&self, closure: f64) -> Result<f64, ModalCouplingError> {
        let mut gradient = [0.0; 2];
        self.storage.gradient(&[-closure, 0.0], &mut gradient);
        let force = finite(-gradient[0])?;
        if force < 0.0 { return Err(invalid("static contact reaction cannot be attractive")); }
        Ok(force)
    }
}

pub(super) fn initialize(
    network: &mut CoupledModalSystem,
    external: &[f64],
    contacts: &[(ModalContact, ModalContactConfig)],
    config: MultiContactConfig,
    gate: &CancelGate,
) -> Result<f64, ModalCouplingError> {
    poll(Some(gate))?;
    network.require_static_initialization(external)?;
    admit_setup(network, contacts.len(), config)?;
    let mut points = Vec::with_capacity(contacts.len());
    for (contact, c) in contacts {
        poll(Some(gate))?;
        let column = contact_column(network, contact, *c)?;
        let storage = ContactStorage::new(Box::new(ZeroStorage), 1, vec![contact.law.clone()])
            .map_err(ModalCouplingError::ContactLaw)?;
        points.push(Point { column, storage, config: *c, gap: contact.law.gaps()[0] });
    }
    let response = StaticResponse::new(network, gate)?;
    // Contact can restrain this prediction; apply physical caps only to the
    // final contacted equilibrium, not a fictional unrestrained load state.
    let free_q = response.solve(network, external, true, false, gate)?;
    let free_x: Vec<f64> = points.iter().map(|p| dot(&p.column, &free_q)).collect::<Result<_,_>>()?;
    let p = points.len();
    let mut compliance = vec![0.0; p*p];
    for j in 0..p {
        poll(Some(gate))?;
        // Rest offsets belong in the free equilibrium, never in a derivative.
        let dq = response.solve(network, &points[j].column, false, false, gate)?;
        for i in 0..p { compliance[i*p+j] = dot(&points[i].column, &dq)?; }
        if compliance[j*p+j] <= 0.0 {
            return Err(invalid("contact attachment needs positive representable static compliance"));
        }
    }
    let mut reactions = vec![0.0; p];
    let sweeps = solve_joint(&compliance, &free_x, |i| points[i].config,
        |i,x| points[i].force(x), config.max_sweeps, &mut reactions, Some(gate))?;
    let mut forces = external.to_vec();
    for (point, reaction) in points.iter().zip(&reactions) {
        if *reaction != 0.0 {
            for (f,b) in forces.iter_mut().zip(&point.column) { *f = finite(*f-b*reaction)?; }
        }
    }
    let q = response.solve(network, &forces, true, true, gate)?;
    let mut energy = network.stage_static_equilibrium(&q, &forces, gate)?;
    // Recompute contact law on ACTUAL rounded states. In addition to the
    // bilateral force-balance gate above, report no equilibrium if one of
    // these contact equations fails. Neither candidate states nor reactions
    // have been published. There is no relaxation-time or dissipation term
    // in this stationary solve, and no discarded startup vibration.
    let mut actual_forces = external.to_vec();
    let mut actual_scale: Vec<f64> = external.iter().map(|f| f.abs()).collect();
    let mut contact_allowance = vec![0.0; network.mode_count()];
    for (point, reaction) in points.iter().zip(&reactions) {
        poll(Some(gate))?;
        let x = extension(&network.candidates, &point.column, 0.0)?;
        let expected = point.force(x)?;
        let residual = finite(reaction-expected)?;
        let tolerance = force_tolerance(*reaction, expected, point.config)?;
        if residual.abs() > tolerance {
            return Err(ModalCouplingError::ContactSolve { residual_n: residual, tolerance_n: tolerance, iterations: sweeps });
        }
        limit("static contact force", expected, point.config.maximum_force_n)?;
        let penetration = finite(x-point.gap)?.max(0.0);
        limit("static contact penetration", penetration, point.config.maximum_penetration_m)?;
        let stored = finite(point.storage.hamiltonian(&[-x, 0.0]))?;
        if stored < 0.0 { return Err(invalid("static contact potential is negative")); }
        energy = finite(energy+stored)?;
        for j in 0..actual_forces.len() {
            let force = finite(point.column[j]*expected)?;
            actual_forces[j] = finite(actual_forces[j]-force)?;
            actual_scale[j] = finite(actual_scale[j]+force.abs())?;
            contact_allowance[j] = finite(contact_allowance[j]+point.column[j].abs()*tolerance)?;
        }
    }
    for (column, link) in network.columns.iter().zip(&network.connections) {
        poll(Some(gate))?;
        let reaction = finite(-link.stiffness_n_m*extension(&network.candidates, column, link.rest_extension_m)?)?;
        for j in 0..actual_forces.len() {
            let force = finite(column[j]*reaction)?;
            actual_forces[j] = finite(actual_forces[j]+force)?;
            actual_scale[j] = finite(actual_scale[j]+force.abs())?;
        }
    }
    let mut j = 0;
    for model in &network.candidates {
        poll(Some(gate))?;
        for (mode, state) in model.modes().iter().zip(model.states()) {
            let actual = finite(mode.angular_frequency_rad_s*mode.angular_frequency_rad_s*state.displacement_m_sqrt_kg)?;
            let scale = finite(actual_scale[j]+actual.abs())?;
            // Contact force tolerances have units N and map through the
            // attachment shapes to N/sqrt(kg). Do not compare these with
            // a dimensionless matrix residual or silently loosen either.
            let tolerance = finite(network.config.solve_relative_tolerance*scale+contact_allowance[j])?;
            if finite(actual-actual_forces[j])?.abs() > tolerance {
                return Err(invalid("actual static contact/network force balance exceeds the mapped force tolerances"));
            }
            j += 1;
        }
    }
    limit("static total energy including contacts", energy, network.config.maximum_total_energy_j)?;
    poll(Some(gate))?;
    std::mem::swap(&mut network.models, &mut network.candidates);
    Ok(energy)
}
