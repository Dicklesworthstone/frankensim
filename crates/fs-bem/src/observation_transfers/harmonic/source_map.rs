//! Cold composition of an original physical field with upstream source laws.
//! The composed field uses the unchanged coherent/statistical, energy,
//! Hermitian pullback and matrix-free Hessian-vector owners.

use super::*;
use fs_math::dd::Dd;

/// Whole-call admission for a new source basis on an already solved field.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicSourceMapLimits {
    pub maximum_sources: usize,
    /// Four complex products per old-source/new-source/receiver tuple.
    /// No sparsity or silent-column discount is used for admission.
    pub maximum_products: usize,
    /// Logical cold peak: supplied map and source identities, new field,
    /// cloned identities/receivers, name-sort indices and identity bytes
    /// hashed. The borrowed original transfer bank and allocator overhead
    /// are excluded. Also bounded by 256 MiB.
    pub maximum_storage_bytes: usize,
}

impl HarmonicSourceMapLimits {
    /// Hard bound shared by adapters that prepare source identities before
    /// invoking the complete composition admission.
    pub const MAXIMUM_STORAGE_BYTES: usize = 256 * 1024 * 1024;
}

impl PreparedHarmonicObservation {
    /// Compose this solved physical field with a complete linear source law.
    ///
    /// `map[new_source * self.sources().len() + old_source]` maps one new
    /// peak coordinate to the indicated original coordinate at this exact
    /// frequency. Its units are old SI unit / new SI unit. Entries retain
    /// complex phase and sign; no conjugation, normalization, passivity
    /// repair or additional physical load is applied. The caller owns the
    /// physical meaning of this upstream law, including impedance/loading.
    ///
    /// All p/u columns are composed in DD before one final round-down. The
    /// result retains the original frequency, medium, receiver geometry and
    /// storage-law capability. Its producer hash binds the original producer,
    /// both ordered physical source identities and the exact complex map.
    /// Preparation is bounded and cancellable; no partial field is published.
    /// Repeated evaluations use the unchanged allocation-free field owners.
    pub fn compose_sources(
        &self,
        sources: &[SourceId],
        map: &[C64],
        limits: HarmonicSourceMapLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        poll(&mut cancelled)?;
        let (n, old, r) = (sources.len(), self.sources.len(), self.receivers.len());
        cap(limits.maximum_sources, MAX_OBSERVATION_INPUTS)?;
        cap(n, limits.maximum_sources)?;
        cap(limits.maximum_products, MAX_PRODUCTS)?;
        if n == 0 || limits.maximum_products == 0 || limits.maximum_storage_bytes == 0 {
            return Err(bad("empty harmonic source-map basis or allowance"));
        }
        let matrix = mul(n, old)?;
        let fields = mul(n, r)?;
        let entries = mul(4, fields)?;
        cap(entries, MAX_OBSERVATION_ENTRIES)?;
        cap(mul(entries, old)?, limits.maximum_products)?;
        if map.len() != matrix {
            return Err(bad(
                "harmonic source-map shape differs from its original basis",
            ));
        }
        // Count complete processing/storage before examining matrix values or
        // cloning names. This also bounds identity hashing on silent maps.
        let mut bytes = mul(matrix, size_of::<C64>())?;
        for (count, size) in [
            (fields, size_of::<[C64; 4]>()),
            (mul(2, n)?, size_of::<SourceId>()),
            (n, size_of::<usize>()),
            (r, size_of::<ReceiverSpec>()),
        ] {
            bytes = add(bytes, mul(count, size)?)?;
        }
        for source in sources {
            poll(&mut cancelled)?;
            bytes = add(bytes, mul(2, source.name().len())?)?;
        }
        for source in &self.sources {
            bytes = add(bytes, source.name().len())?;
        }
        for receiver in &self.receivers {
            bytes = add(bytes, receiver.id.len())?;
        }
        cap(bytes, HarmonicSourceMapLimits::MAXIMUM_STORAGE_BYTES)?;
        cap(bytes, limits.maximum_storage_bytes)?;
        let mut order = allocated(n)?;
        for (index, source) in sources.iter().enumerate() {
            poll(&mut cancelled)?;
            if source.name().is_empty() {
                return Err(bad("harmonic source-map names must be nonempty"));
            }
            order.push(index);
        }
        order.sort_unstable_by(|&a, &b| sources[a].name().cmp(sources[b].name()));
        if order
            .windows(2)
            .any(|pair| sources[pair[0]].name() == sources[pair[1]].name())
        {
            return Err(bad("harmonic source-map names must be distinct"));
        }
        for block in map.chunks(256) {
            poll(&mut cancelled)?;
            if block.iter().any(|&value| !finite(value)) {
                return Err(bad("harmonic source-map coefficients must be finite"));
            }
        }
        let producer_hash = self.source_map_hash(sources, map, &mut cancelled)?;
        let mut transfers = allocated(fields)?;
        for row in map.chunks_exact(old) {
            for receiver in 0..r {
                poll(&mut cancelled)?;
                let mut sum = [[Dd::ZERO; 2]; 4];
                for (old_source, &coefficient) in row.iter().enumerate() {
                    if old_source % 256 == 0 {
                        poll(&mut cancelled)?;
                    }
                    let [ar, ai] = [Dd::from_f64(coefficient.re), Dd::from_f64(coefficient.im)];
                    for (total, &transfer) in sum
                        .iter_mut()
                        .zip(&self.transfers[old_source * r + receiver])
                    {
                        let [br, bi] = [Dd::from_f64(transfer.re), Dd::from_f64(transfer.im)];
                        total[0] = total[0] + ar * br - ai * bi;
                        total[1] = total[1] + ar * bi + ai * br;
                    }
                }
                if sum
                    .iter()
                    .flatten()
                    .any(|value| !value.hi.is_finite() || !value.lo.is_finite())
                {
                    return Err(bad("harmonic source-map composition overflow"));
                }
                let field = sum.map(|parts| C64::new(parts[0].to_f64(), parts[1].to_f64()));
                if field.iter().any(|&z| !finite(z)) {
                    return Err(bad("harmonic source-map composition overflow"));
                }
                transfers.push(field);
            }
        }
        let mut copied_sources = allocated(n)?;
        for source in sources {
            poll(&mut cancelled)?;
            let name = copied_name(source.name())?;
            copied_sources.push(match source {
                SourceId::StructuralMode(_) => SourceId::StructuralMode(name),
                SourceId::ApertureVolumeFlow(_) => SourceId::ApertureVolumeFlow(name),
                SourceId::PointVolumeFlow(_) => SourceId::PointVolumeFlow(name),
                SourceId::PointForce(_) => SourceId::PointForce(name),
                SourceId::OpenCircuitPressure(_) => SourceId::OpenCircuitPressure(name),
                SourceId::PrescribedVolumeFlow(_) => SourceId::PrescribedVolumeFlow(name),
            });
        }
        let mut receivers = allocated(r)?;
        for receiver in &self.receivers {
            poll(&mut cancelled)?;
            receivers.push(ReceiverSpec {
                id: copied_name(&receiver.id)?,
                position_m: receiver.position_m,
            });
        }
        poll(&mut cancelled)?;
        Ok(Self {
            omega: self.omega,
            medium: self.medium,
            sources: copied_sources,
            receivers,
            producer_hash,
            transfers,
            radiation_momentum: self.radiation_momentum,
        })
    }

    fn source_map_hash(
        &self,
        sources: &[SourceId],
        map: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ContentHash, TransferError> {
        let mut hash = DomainHasher::new("fs-bem:harmonic-source-map:v1");
        hash.update(self.producer_hash.as_bytes());
        for value in [self.omega, self.medium.density, self.medium.sound_speed] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        hash.update(&[u8::from(self.radiation_momentum)]);
        for identities in [&self.sources[..], sources] {
            hash.update(&(identities.len() as u64).to_le_bytes());
            for source in identities {
                poll(cancelled)?;
                hash.update(&[source.hash_kind()]);
                hash_name(&mut hash, source.name(), cancelled)?;
            }
        }
        hash.update(&(self.receivers.len() as u64).to_le_bytes());
        for receiver in &self.receivers {
            hash_name(&mut hash, &receiver.id, cancelled)?;
            for value in receiver.position_m {
                hash.update(&value.to_bits().to_le_bytes());
            }
        }
        for block in map.chunks(256) {
            poll(cancelled)?;
            for value in block {
                hash.update(&value.re.to_bits().to_le_bytes());
                hash.update(&value.im.to_bits().to_le_bytes());
            }
        }
        Ok(hash.finalize())
    }
}

fn hash_name(
    hash: &mut DomainHasher,
    name: &str,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), TransferError> {
    hash.update(&(name.len() as u64).to_le_bytes());
    for block in name.as_bytes().chunks(4096) {
        poll(cancelled)?;
        hash.update(block);
    }
    Ok(())
}
fn copied_name(name: &str) -> Result<String, TransferError> {
    let mut copy = String::new();
    copy.try_reserve_exact(name.len())
        .map_err(|_| bad("harmonic source-map name allocation"))?;
    copy.push_str(name);
    Ok(copy)
}
fn allocated<T>(len: usize) -> Result<Vec<T>, TransferError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| bad("harmonic source-map allocation"))?;
    Ok(values)
}
fn add(a: usize, b: usize) -> Result<usize, TransferError> {
    a.checked_add(b)
        .ok_or(bad("harmonic source-map extent overflow"))
}
fn mul(a: usize, b: usize) -> Result<usize, TransferError> {
    a.checked_mul(b)
        .ok_or(bad("harmonic source-map extent overflow"))
}
fn bad(message: &'static str) -> TransferError {
    TransferError::Invalid(message)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), TransferError> {
    if cancelled() {
        Err(TransferError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
