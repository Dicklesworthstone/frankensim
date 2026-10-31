//! Exact bounded bytes for caller-owned checkpoint storage. Decoding is not
//! authentication or physical admission: `restore` must still check fresh
//! endpoints, and a durable consumer must bind its source and original budgets.
use super::*;
use fs_ascent::projected_al::ProjectedAlCheckpoint;

const MAGIC: &[u8; 8] = b"FSMMLCP\x03";

/// Caps apply to the entire payload and cumulative vector entries, including
/// nested lists. They bound decoded storage, not process-wide peak memory.
#[derive(Debug, Clone, Copy)]
pub struct MinimumMaterialCodecBudget3 {
    pub max_bytes: usize,
    pub max_items: usize,
}
impl Default for MinimumMaterialCodecBudget3 {
    fn default() -> Self { Self { max_bytes: 64 * 1024 * 1024, max_items: 4_000_000 } }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MinimumMaterialCodecError3 {
    Invalid(&'static str),
    SizeLimit,
    Allocation,
}
impl std::fmt::Display for MinimumMaterialCodecError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "minimum-material checkpoint bytes refused: {self:?}")
    }
}
impl std::error::Error for MinimumMaterialCodecError3 {}
type Result<T> = std::result::Result<T, MinimumMaterialCodecError3>;
use MinimumMaterialCodecError3::{Allocation, Invalid, SizeLimit};
struct Writer { bytes: Vec<u8>, remaining: usize, items: usize }
impl Writer {
    fn put(&mut self, bytes: &[u8]) -> Result<()> {
        if bytes.len() > self.remaining { return Err(SizeLimit); }
        self.bytes.try_reserve(bytes.len()).map_err(|_| Allocation)?;
        self.bytes.extend_from_slice(bytes); self.remaining -= bytes.len(); Ok(())
    }
}
struct Reader<'a> { bytes: &'a [u8], items: usize }
impl Reader<'_> {
    fn take<const N: usize>(&mut self) -> Result<[u8; N]> {
        let (head, tail) = self.bytes.split_at_checked(N).ok_or(Invalid("truncated checkpoint"))?;
        self.bytes = tail;
        Ok(head.try_into().expect("exact checked byte count"))
    }
}
// Private format implementation: no extensible serialization framework or
// dependency. Every field's order and tag is fixed by this version's layout.
trait Wire: Sized {
    const MIN: usize;
    fn write(&self, w: &mut Writer) -> Result<()>;
    fn read(r: &mut Reader<'_>) -> Result<Self>;
}
impl Wire for fs_blake3::ContentHash {
    const MIN: usize = 32;
    fn write(&self, w: &mut Writer) -> Result<()> { w.put(self.as_bytes()) }
    fn read(r: &mut Reader<'_>) -> Result<Self> { Ok(fs_blake3::ContentHash(r.take()?)) }
}
impl Wire for usize {
    const MIN: usize = 8;
    fn write(&self, w: &mut Writer) -> Result<()> {
        w.put(&u64::try_from(*self).map_err(|_| SizeLimit)?.to_le_bytes())
    }
    fn read(r: &mut Reader<'_>) -> Result<Self> {
        usize::try_from(u64::from_le_bytes(r.take()?)).map_err(|_| SizeLimit)
    }
}
impl Wire for f64 {
    const MIN: usize = 8;
    fn write(&self, w: &mut Writer) -> Result<()> {
        if !self.is_finite() { return Err(Invalid("nonfinite checkpoint number")); }
        w.put(&self.to_bits().to_le_bytes())
    }
    fn read(r: &mut Reader<'_>) -> Result<Self> {
        let v = f64::from_bits(u64::from_le_bytes(r.take()?));
        if v.is_finite() { Ok(v) } else { Err(Invalid("nonfinite checkpoint number")) }
    }
}
impl Wire for bool {
    const MIN: usize = 1;
    fn write(&self, w: &mut Writer) -> Result<()> { w.put(&[u8::from(*self)]) }
    fn read(r: &mut Reader<'_>) -> Result<Self> {
        match r.take()? { [0] => Ok(false), [1] => Ok(true), _ => Err(Invalid("unknown boolean tag")) }
    }
}
impl<T: Wire> Wire for Vec<T> {
    const MIN: usize = 8;
    fn write(&self, w: &mut Writer) -> Result<()> {
        if self.len() > w.items { return Err(SizeLimit); }
        w.items -= self.len(); self.len().write(w)?;
        for item in self { item.write(w)?; } Ok(())
    }
    fn read(r: &mut Reader<'_>) -> Result<Self> {
        let n = usize::read(r)?;
        if n > r.items || n > r.bytes.len() / T::MIN { return Err(SizeLimit); }
        r.items -= n;
        let mut result = Vec::new(); result.try_reserve_exact(n).map_err(|_| Allocation)?;
        for _ in 0..n { result.push(T::read(r)?); } Ok(result)
    }
}
impl<T: Wire> Wire for Option<T> {
    const MIN: usize = 1;
    fn write(&self, w: &mut Writer) -> Result<()> {
        self.is_some().write(w)?;
        match self { Some(v) => v.write(w), None => Ok(()) }
    }
    fn read(r: &mut Reader<'_>) -> Result<Self> {
        if bool::read(r)? { Ok(Some(T::read(r)?)) } else { Ok(None) }
    }
}
impl Wire for ResponseLimitSense3 {
    const MIN: usize = 1;
    fn write(&self, w: &mut Writer) -> Result<()> { matches!(self, Self::AtLeast).write(w) }
    fn read(r: &mut Reader<'_>) -> Result<Self> {
        Ok(if bool::read(r)? { Self::AtLeast } else { Self::AtMost })
    }
}
impl Wire for MinimumMaterialObjective3 {
    const MIN: usize = 1;
    fn write(&self, w: &mut Writer) -> Result<()> {
        match self {
            Self::VolumeFraction => w.put(&[0]),
            Self::MassKg { scale_kg } => { w.put(&[1])?; scale_kg.write(w) },
        }
    }
    fn read(r: &mut Reader<'_>) -> Result<Self> {
        match r.take()? {
            [0] => Ok(Self::VolumeFraction),
            [1] => Ok(Self::MassKg { scale_kg: f64::read(r)? }),
            _ => Err(Invalid("unknown minimum-material objective tag")),
        }
    }
}
macro_rules! fields {
    ($t:ty { $($field:ident : $kind:ty),+ $(,)? }) => {
        impl Wire for $t {
            const MIN: usize = 0 $(+ <$kind as Wire>::MIN)+;
            fn write(&self, w: &mut Writer) -> Result<()> { $(self.$field.write(w)?;)+ Ok(()) }
            fn read(r: &mut Reader<'_>) -> Result<Self> { Ok(Self { $($field: <$kind as Wire>::read(r)?),+ }) }
        }
    };
}
fields!(ProjectedAlWork { iterations: usize, evaluations: usize, multiplier_updates: usize, rejected_trials: usize });
fields!(ProjectedAlOptions { max_dimension: usize, max_evaluations: usize, max_multiplier_updates: usize,
    max_backtracks: usize, tolerance: f64, initial_penalty: f64, max_penalty: f64, inner_tolerance: f64,
    min_spectral_step: f64, max_spectral_step: f64 });
fields!(ProjectedAlSample { objective: f64, gradient: Vec<f64>, constraint: f64, constraint_gradient: Vec<f64> });
fields!(ProjectedAlCheckpoint { point: Vec<f64>, sample: ProjectedAlSample, work: ProjectedAlWork,
    multiplier: f64, penalty: f64, spectral_step: f64, inner_tolerance: f64, previous_violation: Option<f64> });
fields!(ResponseLimitOptions3 { sharpness: f64, max_cases: usize, max_limits: usize });
fields!(MinimumMaterialOptions3 { limits: ResponseLimitOptions3, density_floor: f64, volume_scale: f64,
    objective: MinimumMaterialObjective3, optimizer: ProjectedAlOptions });
fields!(MinimumMaterialLimit3 { reaction: bool, sense: ResponseLimitSense3, bound: f64, scale: f64 });
fields!(MinimumMaterialEndpoint3 { rho: Vec<f64>, projected_rho: Vec<f64>, scales: Vec<f64>,
    volume_fraction: f64, mass_kg: Option<f64>, objective: f64, constraint: f64, worst_residual: f64,
    measurements: Vec<Vec<f64>>, residuals: Vec<Vec<f64>> });
fields!(MinimumMaterialIteration3 { iteration: usize, volume_fraction: f64, mass_kg: Option<f64>,
    objective: f64, constraint: f64, worst_residual: f64 });
fields!(MinimumMaterialCheckpoint3 { optimizer: ProjectedAlCheckpoint, options: MinimumMaterialOptions3,
    material_body_load_fingerprint: Option<fs_blake3::ContentHash>,
    limits: Vec<Vec<MinimumMaterialLimit3>>, prescribed: Vec<bool>, reference_cell_masses_kg: Option<Vec<f64>>,
    accepted: MinimumMaterialEndpoint3,
    best_feasible: Option<MinimumMaterialEndpoint3>, history: Vec<MinimumMaterialIteration3>, restoration_evaluations: usize });

impl MinimumMaterialCheckpoint3 {
    /// Encode finite values with exact IEEE bits and little-endian u64 counts.
    /// No authentication or source/budget authority is conferred by these bytes.
    pub fn to_bytes(&self, budget: MinimumMaterialCodecBudget3) -> Result<Vec<u8>> {
        let mut w = Writer { bytes: Vec::new(), remaining: budget.max_bytes, items: budget.max_items };
        w.put(MAGIC)?; self.write(&mut w)?; Ok(w.bytes)
    }
    /// Reject unknown versions/tags, nonfinite values, oversized lengths,
    /// truncated payloads and trailing data before returning a decoded value.
    /// The returned data must still pass `MinimumMaterialStudy3::restore`.
    pub fn from_bytes(bytes: &[u8], budget: MinimumMaterialCodecBudget3) -> Result<Self> {
        if bytes.len() > budget.max_bytes { return Err(SizeLimit); }
        let mut r = Reader { bytes, items: budget.max_items };
        if &r.take::<8>()? != MAGIC { return Err(Invalid("unknown checkpoint version")); }
        let result = Self::read(&mut r)?;
        if !r.bytes.is_empty() { return Err(Invalid("trailing checkpoint data")); }
        Ok(result)
    }
}
