// Shared cold source ingress for live instruments. Only frame unchanged
// canonical bytes; species matching, thermodynamics, transport mixing and
// source validity remain in fs-material's read_source_packet.
export function gasSourcePacket(record, sources) {
  if (sources === undefined) return undefined;
  if (typeof record !== 'string' || !/^\s*thermochemical\s/.test(record)
      || !sources || !Array.isArray(sources.components)
      || sources.components.length < 1 || sources.components.length > 128
      || !['frozen', 'liquid-water'].includes(sources.phase)
      || !(sources.waterSpecies === null || (typeof sources.waterSpecies === 'string'
        && sources.waterSpecies.trim() === sources.waterSpecies && sources.waterSpecies.length > 0))
      || (sources.phase === 'liquid-water' && sources.waterSpecies === null)) {
    throw new Error('Declare a thermochemical gas record, 1..128 source triples, explicit water identity (or null), and phase policy.');
  }
  const parts = [], lengths = [], limit = 16 * 1024 * 1024;
  let total = 0;
  for (const component of sources.components) {
    for (const name of ['speciesPack', 'nasa9Pack', 'transportPack']) {
      const bytes = component?.[name];
      if (!(bytes instanceof Uint8Array) || !bytes.length || bytes.length > limit) {
        throw new Error('Each gas source triple needs three nonempty canonical packs, at most 16 MiB each.');
      }
      total += bytes.length;
      if (total > 4 * limit) throw new Error('Gas source packs exceed the aggregate 64 MiB budget.');
      parts.push(bytes); lengths.push(bytes.length);
    }
  }
  const bytes = new Uint8Array(total);
  let offset = 0;
  for (const part of parts) { bytes.set(part, offset); offset += part.length; }
  return [bytes, new Uint32Array(lengths), sources.waterSpecies ?? '', sources.phase];
}
