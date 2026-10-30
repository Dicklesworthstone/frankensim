//! The selected shell receiver on the actual connected harmonic response.
//! Raw exp(-i omega t) producer fields retain physical propagation phase. FIR
//! latency, taper and sample history do not enter this continuous observation.
//! Geometry, source traces and mouth inputs are the playback owners' originals;
//! no additional air reaction or loss is debited by observing them.

use super::{ShellObserverError, ShellObserverOptions, invalid, loaded_observation_inputs, poll};
use super::one_way::{self, OneWayShellObserverOptions, ShellAcousticSkin, skin_identity, skin_traces};
use crate::bowed_string::runtime::reciprocal::bank::harmonic::ShellBridgeHarmonicResponse;
use crate::bowed_string::runtime::reciprocal::bank::multiple::MultiBowStringBank;
use fs_bem::helmholtz::{self, Medium};
use fs_bem::near_field::Options as ReceiverQuadrature;
use fs_bem::observation_transfers::{
    BemObservation, ReciprocalBemObservation, ReceiverSpec, SourceId, TransferError,
};
use fs_bem::panel3d::SpherePanels;
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_math::c64::C64;

/// Reuse one explicitly selected playback source/receiver declaration. Only
/// its physical source, pass band and numerical boundary controls are used.
pub enum ShellHarmonicReceiverSelection<'a> {
    Loaded { options: &'a ShellObserverOptions, quadrature: ReceiverQuadrature },
    Skin { options: &'a OneWayShellObserverOptions, skin: &'a ShellAcousticSkin },
}

/// Physical scope of the raw producer. The skin observer retains its original
/// ruled collocation formulation; it is not relabeled as its reciprocal load.
#[derive(Clone, Debug)]
pub struct ShellHarmonicReceiverReport {
    pub receiver: ReceiverSpec,
    pub gas: GasState,
    pub formulation: &'static str,
    pub operator_identity: ContentHash,
    pub loading_identity: Option<ContentHash>,
    pub skin_identity: Option<ContentHash>,
    pub body_modes: usize,
    pub opening_ports: usize,
    pub surface_panels: usize,
    /// Structural modes receiving the actually retained exterior reaction.
    pub loaded_body_modes: usize,
    /// Other structural modes are observed but have no exterior feedback.
    /// This count does not establish that omitted reaction is negligible.
    pub one_way_body_modes: usize,
    pub pass_band_upper_hz: f64,
    /// Conservative source solve and receiver quadrature reservation per query.
    /// Coherent projection of each returned drive is reserved separately.
    pub observation_work_per_frequency: u64,
}

/// One force experiment, with the selected receiver separate from the native
/// compact/modal pressure retained inside `body.released`.
#[derive(Clone, Debug)]
pub struct ShellReceiverHarmonicResponse {
    pub body: ShellBridgeHarmonicResponse,
    pub receiver_pressure_pa: C64,
    pub receiver_gradient_pa_per_m: [C64; 3],
}

enum Kind {
    Loaded(ShellObserverOptions, ReceiverQuadrature),
    Skin(OneWayShellObserverOptions),
}

/// Geometry and complete physical input map of one immutable borrowed bank.
/// Preparing the receiver borrows this geometry, so no self-reference, copied
/// mechanics or independently supplied same-sized response can be cross-wired.
pub struct ShellHarmonicReceiver<'a> {
    bank: &'a MultiBowStringBank,
    surface: SpherePanels,
    fields: Vec<Vec<f64>>,
    sources: Vec<SourceId>,
    coordinates: Vec<(usize, f64)>,
    bridge_count: usize,
    frequency_samples: usize,
    kind: Kind,
    report: ShellHarmonicReceiverReport,
}

impl<'a> ShellHarmonicReceiver<'a> {
    /// Admit the same shell or skin and every original modal/mouth source.
    /// `frequency_samples` bounds the declared observation sweep. Each query
    /// additionally charges the caller's combined mechanics/observation work.
    /// The skin's modal omitted-load adequacy sweep is not repeated here; the
    /// report retains actual loading scope without promoting one-way accuracy.
    pub fn new(bank: &'a MultiBowStringBank, selection: ShellHarmonicReceiverSelection<'_>,
        frequency_samples: usize, gate: &CancelGate) -> Result<Self, ShellObserverError> {
        poll(gate)?;
        if frequency_samples == 0 { return Err(invalid("harmonic receiver needs a nonempty sweep")); }
        let source = bank.source().bank();
        let components = bank.system().components();
        let body_component = source.strings().len();
        let body_modes = components.get(body_component)
            .ok_or_else(||invalid("harmonic receiver lost its original body"))?.modes().len();
        let mut parts = Vec::new();
        let bridge_count = if let Some(image) = source.material_shell_image() {
            parts.push((image.source().mesh(),image.physical_modes()));
            image.force_ports().len()
        } else if let Some(image) = source.material_shell_assembly_image() {
            for (part, rows) in image.source().components.iter().zip(image.physical_modes()) {
                parts.push((part.specimen.mesh(),rows.as_slice()));
            }
            image.force_ports().len()
        } else { return Err(invalid("harmonic receiver requires its original material shell")); };
        if body_modes == 0 || bridge_count == 0
            || parts.iter().any(|(_, rows)|rows.len()!=body_modes) {
            return Err(invalid("harmonic receiver must retain every original shell mode"));
        }
        let gas = *source.gas();
        gas.validate_coherence().map_err(|_|invalid("harmonic receiver gas is incoherent"))?;
        let loading = source.body_radiation().map(|binding| {
            let report = binding.material_shell_bem()
                .ok_or_else(||invalid("raw or planar air loads have no shell receiver source binding"))?;
            if binding.gas()!=source.gas() || report.gas!=gas
                || binding.basis_identity()!=report.basis_identity {
                return Err(invalid("harmonic receiver load lost its original basis or gas"));
            }
            Ok((binding,report))
        }).transpose()?;
        let (surface, fields, sources, addresses, kind, receiver, pass_band_upper_hz,
            formulation, skin_hash, loaded_body_modes) = match selection {
            ShellHarmonicReceiverSelection::Loaded {options,quadrature} => {
                options.validate(source.options().sample_rate_hz)?;
                quadrature.validate().map_err(TransferError::Helmholtz)?;
                let (binding, loading) = loading
                    .ok_or_else(||invalid("loaded harmonic receiver requires its original BEM load"))?;
                if source.material_cavity().is_some() || loading.skin_identity.is_some()
                    || loading.body_modes!=body_modes || loading.one_way_body_modes!=0
                    || options.pass_band_upper_hz>loading.options.maximum_hz {
                    return Err(invalid("loaded harmonic receiver requires the complete matching shell source and band"));
                }
                let expected = source.radiation_body_coordinates(body_modes)
                    .map_err(ShellObserverError::Dynamics)?;
                if expected!=binding.source_body_coordinates()
                    || expected.len()!=body_modes+loading.opening_ports {
                    return Err(invalid("harmonic receiver lost its complete loaded body/mouth map"));
                }
                let (surface, mut fields, _) = source.shell_surface_traces(body_modes,
                    loading.options.maximum_surface_panels,options.maximum_sweep_terms,
                    frequency_samples,gate).map_err(ShellObserverError::Shell)?;
                if surface.areas().len()!=loading.surface_panels {
                    return Err(invalid("harmonic receiver geometry differs from its loaded shell"));
                }
                let (sources, addresses) = loaded_observation_inputs(source,body_modes,
                    loading,&expected,&mut fields)?;
                (surface,fields,sources,addresses,Kind::Loaded(options.clone(),quadrature),
                    options.receiver.clone(),options.pass_band_upper_hz,
                    "reciprocal_p1_p0",None,body_modes)
            }
            ShellHarmonicReceiverSelection::Skin {options,skin} => {
                options.validate(source.options().sample_rate_hz)?;
                if source.spatial_cavity().is_some() {
                    return Err(invalid("skin harmonic receiver has no spatial cavity mouth binding"));
                }
                let identity = skin_identity(skin);
                let loaded = if let Some((binding, loading)) = loading {
                    let neck = usize::from(source.material_cavity().is_some());
                    let expected: Vec<_> = (0..loading.body_modes)
                        .chain((0..neck).map(|_|body_modes)).collect();
                    if loading.skin_identity!=Some(identity)
                        || loading.body_modes+loading.one_way_body_modes!=body_modes
                        || loading.opening_ports!=neck || expected!=binding.source_body_coordinates() {
                        return Err(invalid("skin harmonic receiver differs from its original loaded skin or mode prefix"));
                    }
                    loading.body_modes
                } else { 0 };
                let cavity = source.material_cavity();
                if cavity.is_some()!=!skin.apertures.is_empty()
                    || (loaded==0 && cavity.is_some_and(|c|
                        c.cavity.loss()==crate::cavity_modal::CavityNeckLoss::ExteriorMouth)) {
                    return Err(invalid("skin harmonic receiver needs the original physical mouth and exterior load ownership"));
                }
                let traces = skin_traces(&parts,body_modes,skin,options.maximum_skin_panels,
                    options.maximum_anchor_distance_m,gate)?;
                if loading.is_some_and(|(_,l)|l.surface_panels!=traces.surface.areas().len()) {
                    return Err(invalid("harmonic receiver panel count differs from its loaded skin"));
                }
                let (sources,addresses,fields) = one_way::observation_inputs(
                    source,components,body_modes,skin,&traces)?;
                let work = one_way::frequency_terms(traces.surface.areas().len(),sources.len())
                    .and_then(|w|w.checked_mul(frequency_samples as u64))
                    .ok_or_else(||invalid("harmonic skin sweep work overflow"))?;
                if work>options.maximum_boundary_terms {
                    return Err(invalid("harmonic skin sweep exceeds its original boundary allowance"));
                }
                (traces.surface,fields,sources,addresses,Kind::Skin(options.clone()),
                    options.receiver.clone(),options.pass_band_upper_hz,
                    "ruled_burton_miller",Some(identity),loaded)
            }
        };
        let sizes: Vec<_> = components.iter().map(|c|c.modes().len()).collect();
        let coordinates = flatten_coordinates(&sizes,&addresses)?;
        let report = ShellHarmonicReceiverReport {
            receiver,gas,formulation,operator_identity:ContentHash([0;32]),
            loading_identity:loading.map(|(_,l)|l.loading_identity),skin_identity:skin_hash,
            body_modes,opening_ports:sources.len()-body_modes,surface_panels:surface.areas().len(),
            loaded_body_modes,one_way_body_modes:body_modes-loaded_body_modes,
            pass_band_upper_hz,observation_work_per_frequency:0,
        };
        poll(gate)?;
        Ok(Self {bank,surface,fields,sources,coordinates,bridge_count,frequency_samples,kind,report})
    }

    /// Retain one original boundary geometry/trace preparation across the sweep.
    /// The reciprocal source caches static moments; each frequency still uses
    /// its own original boundary operator and one batched source factorization.
    pub fn prepare(&self, gate: &CancelGate)
        -> Result<PreparedShellHarmonicReceiver<'_, 'a>, ShellObserverError> {
        poll(gate)?;
        let gas = self.report.gas;
        let medium = Medium {density:gas.density,sound_speed:gas.sound_speed};
        let mut report = self.report.clone();
        let (owner, work) = match &self.kind {
            Kind::Loaded(options,quadrature) => {
                let loading = self.bank.source().bank().body_radiation()
                    .and_then(|b|b.material_shell_bem())
                    .ok_or_else(||invalid("harmonic receiver lost its immutable BEM binding"))?;
                let mut allowance = loading.options;
                allowance.maximum_sweep_terms=options.maximum_sweep_terms;
                let owner = ReciprocalBemObservation::new(&self.surface,medium,&self.sources,
                    &self.fields,core::slice::from_ref(&report.receiver),options.minimum_clearance_m,
                    core::f64::consts::TAU*report.pass_band_upper_hz,
                    allowance.reciprocal_budget(self.frequency_samples,loading.reciprocal.angular_refinement)
                        .map_err(ShellObserverError::Shell)?,*quadrature,||gate.is_requested())?;
                let work = owner.frequency_work().map_err(TransferError::Helmholtz)?
                    .checked_add(quadrature.maximum_kernel_evaluations as u64)
                    .ok_or_else(||invalid("harmonic reciprocal receiver work overflow"))?;
                report.operator_identity=owner.operator_hash();
                (FrequencyOwner::Loaded(owner),work)
            }
            Kind::Skin(options) => {
                let owner = BemObservation::new_ruled_burton_miller(&self.surface,medium,
                    &self.sources,&self.fields,core::slice::from_ref(&report.receiver),
                    options.minimum_clearance_m,options.panel_rule)?;
                let work = one_way::frequency_terms(self.surface.areas().len(),self.sources.len())
                    .ok_or_else(||invalid("harmonic skin receiver work overflow"))?;
                report.operator_identity=owner.operator_hash();
                (FrequencyOwner::Skin(owner,options.workers),work)
            }
        };
        report.observation_work_per_frequency=work;
        poll(gate)?;
        Ok(PreparedShellHarmonicReceiver {source:self,owner,report})
    }
}

enum FrequencyOwner<'a> {
    Loaded(ReciprocalBemObservation<'a>),
    Skin(BemObservation<'a>,usize),
}

/// One receiver bound to one original bank and source geometry. Queries are
/// read-only; each complete bridge matrix shares one BEM frequency sample.
pub struct PreparedShellHarmonicReceiver<'a,'b> {
    source: &'a ShellHarmonicReceiver<'b>,
    owner: FrequencyOwner<'a>,
    report: ShellHarmonicReceiverReport,
}
impl PreparedShellHarmonicReceiver<'_, '_> {
    pub fn report(&self) -> &ShellHarmonicReceiverReport { &self.report }

    fn mechanics_work(&self,omega:f64,drives:usize,maximum_work:usize)
        ->Result<usize,ShellObserverError> {
        if !omega.is_finite() || omega<=0.0
            || omega>core::f64::consts::TAU*self.report.pass_band_upper_hz {
            return Err(invalid("harmonic receiver frequency exceeds its original physical pass band"));
        }
        remaining_mechanics_work(maximum_work,self.report.observation_work_per_frequency,
            self.source.coordinates.len(),drives)
    }

    fn sample(&self,omega:f64,gate:&CancelGate)->Result<Vec<C64>,ShellObserverError> {
        poll(gate)?;
        let row = match &self.owner {
            FrequencyOwner::Loaded(owner) => owner.frequency_sample(omega,||gate.is_requested()),
            FrequencyOwner::Skin(owner,workers) =>
                helmholtz::with_dense_workers(*workers,||owner.frequency_sample(omega)),
        }.map_err(TransferError::Helmholtz)?;
        poll(gate)?;
        Ok(row)
    }

    fn observe(&self,body:ShellBridgeHarmonicResponse,row:&[C64],gate:&CancelGate)
        ->Result<ShellReceiverHarmonicResponse,ShellObserverError> {
        poll(gate)?;
        let field = project_field(row,&self.source.coordinates,
            &body.released.dynamics.velocity_m_sqrt_kg_per_s)?;
        Ok(ShellReceiverHarmonicResponse {body,receiver_pressure_pa:field[0],
            receiver_gradient_pa_per_m:[field[1],field[2],field[3]]})
    }

    pub fn bridge_harmonic_response(&self,omega:f64,bridge:usize,force:C64,
        maximum_work:usize,gate:&CancelGate)->Result<ShellReceiverHarmonicResponse,ShellObserverError> {
        let work = self.mechanics_work(omega,1,maximum_work)?;
        let body = self.source.bank.material_shell_bridge_harmonic_response(omega,bridge,force,work,gate)
            .map_err(ShellObserverError::Dynamics)?;
        let row = self.sample(omega,gate)?;
        self.observe(body,&row,gate)
    }

    /// Coherent physical force phasors on all declared bridges, with all cross
    /// terms in mechanics, heat and receiver interference retained.
    pub fn bridge_forces_harmonic_response(&self,omega:f64,forces:&[C64],
        maximum_work:usize,gate:&CancelGate)->Result<ShellReceiverHarmonicResponse,ShellObserverError> {
        let work = self.mechanics_work(omega,1,maximum_work)?;
        let body = self.source.bank.material_shell_bridge_forces_harmonic_response(omega,forces,work,gate)
            .map_err(ShellObserverError::Dynamics)?;
        let row = self.sample(omega,gate)?;
        self.observe(body,&row,gate)
    }

    /// Unit physical bridge-force columns in declaration order. The existing
    /// mechanics shares its factorization; this receiver shares its entire
    /// forcing-frequency source solve and field rows across all columns. The
    /// skin producer remains synchronous within each original boundary sample;
    /// reciprocal-source cancellation also reaches its LU/quadrature rows.
    pub fn bridge_harmonic_matrix(&self,omega:f64,maximum_work:usize,gate:&CancelGate)
        ->Result<Vec<ShellReceiverHarmonicResponse>,ShellObserverError> {
        let work = self.mechanics_work(omega,self.source.bridge_count,maximum_work)?;
        let bodies = self.source.bank.material_shell_bridge_harmonic_matrix(omega,work,gate)
            .map_err(ShellObserverError::Dynamics)?;
        let row = self.sample(omega,gate)?;
        bodies.into_iter().map(|body|self.observe(body,&row,gate)).collect()
    }
}

fn flatten_coordinates(sizes:&[usize],addresses:&[(usize,usize,f64)])
    ->Result<Vec<(usize,f64)>,ShellObserverError> {
    let mut offsets=Vec::with_capacity(sizes.len());let mut total=0usize;
    for &size in sizes {
        offsets.push(total);
        total=total.checked_add(size).ok_or_else(||invalid("harmonic source coordinate extent overflow"))?;
    }
    addresses.iter().map(|&(component,coordinate,scale)| {
        if sizes.get(component).is_none_or(|&size|coordinate>=size) || !scale.is_finite() {
            return Err(invalid("harmonic source lost its original component coordinate or physical scale"));
        }
        Ok((offsets[component]+coordinate,scale))
    }).collect()
}

fn remaining_mechanics_work(maximum:usize,observation:u64,sources:usize,drives:usize)
    ->Result<usize,ShellObserverError> {
    // Four complex products/sums plus physical input scaling per source/drive.
    let projection=sources.checked_mul(drives).and_then(|v|v.checked_mul(40));
    let reserved=usize::try_from(observation).ok().and_then(|v|projection.and_then(|p|v.checked_add(p)));
    if sources==0 || drives==0 {return Err(invalid("harmonic receiver needs sources and force columns"));}
    reserved.and_then(|r|maximum.checked_sub(r)).filter(|&v|v!=0)
        .ok_or_else(||invalid("complete harmonic response exceeds its combined mechanics/receiver work allowance"))
}

// Source-major p,dp/dx,dp/dy,dp/dz rows of the original exp(-iwt) producer.
// Complex source velocities are multiplied without conjugation; cancellation
// between structural modes and mouth flows is part of the physical field.
fn project_field(row:&[C64],coordinates:&[(usize,f64)],velocity:&[C64])
    ->Result<[C64;4],ShellObserverError> {
    if row.len()!=coordinates.len().checked_mul(4).unwrap_or(usize::MAX) {
        return Err(invalid("harmonic receiver lost its complete source-major field rows"));
    }
    let mut field=[C64::ZERO;4];
    for (transfer,&(index,scale)) in row.chunks_exact(4).zip(coordinates) {
        let v=velocity.get(index).ok_or_else(||invalid("harmonic receiver lost its original source velocity"))?.scale(scale);
        for (out,&h) in field.iter_mut().zip(transfer) { *out=*out+h*v; }
    }
    if field.iter().any(|z|!z.re.is_finite() || !z.im.is_finite()) {
        return Err(invalid("harmonic receiver pressure or gradient is not representable"));
    }
    Ok(field)
}

#[cfg(test)]
mod tests;
