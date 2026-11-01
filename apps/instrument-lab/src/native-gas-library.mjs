// Original source selection only. Gas laws, source identities and physical
// domains remain in the native material compiler and fs-material registry.
import library from './native-gas-library.json' with {type:'json'};
import {MAX_ATTACHMENTS, MAX_ATTACHMENT_FILE_BYTES, MAX_ATTACHMENT_BYTES, admitAttachmentPath} from './native-experiment.mjs';

export const NATIVE_GAS_LIBRARY = library;

/** Load all six original files before publishing any part of a gas source. */
export async function loadNativeGasSource(id, readSource) {
  const entry = library.find(source => source.id === id);
  if (!entry) throw new Error('Choose a gas with complete species, NASA-9 and transport sources.');
  const prefix = `gas-sources/${entry.id}`;
  const descriptors = [
    ['species', entry.speciesSource, 'species.tsv'],
    ['nasa9', entry.nasa9Source, 'nasa9.tsv'],
    ['transport', entry.transportSource, 'properties.tsv'],
  ];
  const attachments = await Promise.all(descriptors.flatMap(([kind, directory, table]) =>
    ['manifest.tsv', table].map(async name => {
      const text = await readSource(`data/matdb/seed-v1/${directory}/${name}`);
      if (typeof text !== 'string') throw new Error('Original gas source is not UTF-8 text.');
      const bytes = new TextEncoder().encode(text);
      if (!bytes.length || bytes.length > MAX_ATTACHMENT_FILE_BYTES)
        throw new Error('Original gas source exceeds the attachment limit.');
      return {path:`${prefix}/${kind}/${name}`, bytes};
    })));
  return {attachments,
    materialManifests:descriptors.map(([kind]) => ({
      manifestPath:`${prefix}/${kind}/manifest.tsv`, packPath:`${prefix}/${kind}.pack`,
    })),
    component:{speciesPath:`${prefix}/species.pack`, nasa9Path:`${prefix}/nasa9.pack`,
      transportPath:`${prefix}/transport.pack`},
  };
}

/** Add sources without replacing an uploaded or edited file/material binding. */
export function mergeNativeGasSource(attachments, materialManifests, added, maximumAttachments = MAX_ATTACHMENTS) {
  const files = [...attachments], bindings = [...(materialManifests ?? [])];
  for (const candidate of added.attachments) {
    const previous = files.find(file => file.path === candidate.path);
    if (previous) {
      if (previous.bytes.length !== candidate.bytes.length
          || previous.bytes.some((value, index) => value !== candidate.bytes[index]))
        throw new Error(`Gas source ${candidate.path} already contains different data. Rename it or keep your existing source.`);
    } else files.push(candidate);
  }
  for (const candidate of added.materialManifests) {
    const previous = bindings.find(binding => binding.packPath === candidate.packPath);
    if (files.some(file => file.path === candidate.packPath)
        || (previous && (previous.manifestPath !== candidate.manifestPath || previous.materialState !== undefined)))
      throw new Error(`Gas pack ${candidate.packPath} is already supplied by another source. Keep or rename that binding explicitly.`);
    if (!previous) bindings.push(candidate);
  }
  const paths = new Set();
  for (const path of [...files.map(file => file.path), ...bindings.map(binding => binding.packPath)])
    admitAttachmentPath(path, paths);
  // Enforce the existing byte limits even on a host without a native renderer.
  // The count includes the three compiled packs as well as their original files.
  const limit = Math.min(maximumAttachments, MAX_ATTACHMENTS);
  if (!Number.isSafeInteger(limit) || limit < 1 || files.length + bindings.length > limit)
    throw new Error(`Original sources and compiled packs exceed the ${limit} file allowance. Update the local renderer or remove unused sources.`);
  if (files.some(file => !file.bytes.length || file.bytes.length > MAX_ATTACHMENT_FILE_BYTES)
      || files.reduce((sum, file) => sum + file.bytes.length, 0) > MAX_ATTACHMENT_BYTES)
    throw new Error('Gas sources exceed the original attachment byte limits.');
  return {attachments:files, materialManifests:bindings};
}
