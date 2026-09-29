// Nephron segments (paper: "Kidney: Nephrons (Podocytes, Proximal Tubule
// Cells, Distal Tubule Cells, Collecting Duct Cells), Renal Corpuscle") and a
// water-reabsorption profile along the tubule.
//
// Fractions of the filtered water reabsorbed per segment are standard
// physiology (Boron & Boulpaep, Medical Physiology; Guyton & Hall):
//   proximal tubule ~65 %, descending thin limb ~15 %, ascending limb ~0 %
//   (water-impermeable, NaCl pumped out), distal convoluted tubule ~5 %,
//   collecting duct 0 .. ~14 % of the filtrate depending on ADH (aquaporin-2).
// Within a segment the remaining volume decays exponentially.

export const NEPHRON_SEGMENTS = [
  { key: 'pct', name: 'Proximal convoluted tubule', cells: 'Proximal tubule cells (brush border)', reabsorbs: '≈65 % water, Na⁺, glucose, amino acids',
    water: 0.65, color: 0x62d2a2, radius: 9, cellsPerRing: 11, width: 2.9, height: 4.4, basalNucleus: true },
  { key: 'tdl', name: 'Loop of Henle · thin descending limb', cells: 'Squamous thin-limb cells', reabsorbs: '≈15 % water (aquaporin-1)',
    water: 0.15, color: 0x7fb8ff, radius: 6.5, cellsPerRing: 8, width: 2.8, height: 1.6, basalNucleus: false },
  { key: 'tal', name: 'Loop of Henle · thick ascending limb', cells: 'NKCC2-expressing epithelium', reabsorbs: '≈25 % NaCl, no water',
    water: 0.0, color: 0xc38bff, radius: 7.5, cellsPerRing: 10, width: 2.6, height: 3.2, basalNucleus: false },
  { key: 'dct', name: 'Distal convoluted tubule', cells: 'Distal tubule cells (ion balance)', reabsorbs: '≈5 % NaCl / water, Ca²⁺',
    water: 0.05, color: 0xffb454, radius: 7.5, cellsPerRing: 10, width: 2.6, height: 3.4, basalNucleus: false },
  { key: 'cd', name: 'Collecting duct', cells: 'Principal & intercalated cells', reabsorbs: 'water on demand (ADH), acid–base',
    water: 0.14, color: 0xff6fae, radius: 10, cellsPerRing: 12, width: 3.0, height: 3.8, basalNucleus: false, adhRegulated: true },
];

// arc-length boundaries of the schematic centre-line used by the scene
export const SEGMENT_BOUNDS = [0, 7 / 20, 10 / 20, 14 / 20, 18 / 20, 1];

/**
 * Fraction of the filtered volume still in the lumen at position u in [0, 1].
 * @param {number} adh 0 (diabetes insipidus) .. 1 (maximal antidiuresis)
 */
export function filtrateProfile(adh = 0.5, segments = NEPHRON_SEGMENTS, bounds = SEGMENT_BOUNDS) {
  const ends = [];
  let remaining = 1;
  for (const s of segments) {
    const start = remaining;
    if (s.adhRegulated) remaining = start * (1 - 0.95 * adh);
    else remaining = start - s.water;
    ends.push([start, Math.max(remaining, 1e-4)]);
  }
  return (u) => {
    const x = Math.min(Math.max(u, 0), 1);
    let i = bounds.findIndex((b, k) => x >= b && x < bounds[k + 1]);
    if (i < 0) i = segments.length - 1;
    const [a, b] = ends[i];
    const f = (x - bounds[i]) / (bounds[i + 1] - bounds[i]);
    return a * Math.pow(b / a, f);
  };
}

export const REACTIONS = [
  { id: 'kidney.filtration', order: 4, pathway: 'Kidney', reactants: ['Blood plasma'], products: ['Glomerular filtrate'], modifiers: ['Podocyte', 'Glomerular endothelium'], label: 'Glomerular filtration' },
  { id: 'kidney.pct', order: 4, pathway: 'Kidney', reactants: ['Glomerular filtrate'], products: ['Blood plasma'], modifiers: ['Proximal tubule cell'], label: 'Proximal reabsorption' },
  { id: 'kidney.loop', order: 4, pathway: 'Kidney', reactants: ['Tubular NaCl'], products: ['Interstitial NaCl'], modifiers: ['NKCC2'], label: 'Countercurrent NaCl transport' },
  { id: 'kidney.adh', order: 4, pathway: 'Kidney', reactants: ['Tubular water'], products: ['Blood plasma'], modifiers: ['ADH', 'Aquaporin-2'], label: 'ADH-regulated water reabsorption' },
];
