// Scene registry. Scenes load lazily so the first paint only pays for one.
export const SCENES = [
  {
    id: 'colony', group: 'Real 3D microscopy', order: 'Real data',
    title: 'Stem-cell colony · digital twin',
    summary: 'Real 3D confocal stack → segmented cells → growing 3D model',
    load: () => import('./colony.js'),
  },
  {
    id: 'kidney', group: 'Real 3D microscopy', order: 'Real data',
    title: 'Kidney tissue · organ-level volume',
    summary: 'Nephron tubules & glomeruli, 3D confocal',
    load: () => import('./kidney.js'),
  },
  {
    id: 'signaling', group: 'Paper stages · orders 1–4', order: 'Order 1',
    title: 'Stage 1 · Notch, Wnt & Hedgehog',
    summary: 'Receptors, cleavage, nuclear translocation, lateral inhibition',
    load: () => import('./signaling.js'),
  },
  {
    id: 'hematopoiesis', group: 'Paper stages · orders 1–4', order: 'Order 2',
    title: 'Stage 2 · Progenitor commitment',
    summary: 'MAPK/ERK, JAK/STAT, GATA1–PU.1, blood lineages',
    load: () => import('./hematopoiesis.js'),
  },
  {
    id: 'differentiation', group: 'Paper stages · orders 1–4', order: 'Order 3',
    title: 'Stage 3 · Muscle, neurons & feedback',
    summary: 'Myoblast fusion, neurite guidance, hormonal loops',
    load: () => import('./differentiation.js'),
  },
  {
    id: 'morphogenesis', group: 'Paper stages · orders 1–4', order: 'Order 4',
    title: 'Stage 4 · Morphogenesis',
    summary: 'Turing organoids, morphogen gradients, sorting, ECM',
    load: () => import('./morphogenesis.js'),
  },
  {
    id: 'proteins', group: 'Molecular & systems', order: 'Molecular',
    title: 'Protein foundry · hemoglobin & collagen',
    summary: 'Transcription → translation → folding → assembly',
    load: () => import('./proteins.js'),
  },
  {
    id: 'landscape', group: 'Molecular & systems', order: 'Epigenetics',
    title: 'Epigenetic landscape',
    summary: 'Waddington surface from −ln P, DNMT / HAT / HDAC',
    load: () => import('./landscape.js'),
  },
  {
    id: 'network', group: 'Molecular & systems', order: 'Orders 1–4',
    title: 'Reaction network · PCA embedding',
    summary: 'Every modelled reaction, grouped by order',
    load: () => import('./network.js'),
  },
  {
    id: 'atlas', group: 'Molecular & systems', order: 'Hierarchy',
    title: 'Cellular map · cells to body',
    summary: 'Lineages, tissues, organs, organ systems',
    load: () => import('./atlas.js'),
  },
  {
    id: 'gallery', group: 'Reference', order: 'Microscopy',
    title: 'Micrograph gallery · animated variations',
    summary: 'Real images, measurements and their animated models',
    load: () => import('./gallery.js'),
  },
];
