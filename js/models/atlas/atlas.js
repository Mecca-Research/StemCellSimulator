// The paper's "Cellular Map of the Human Body": specialised cell types grouped
// by lineage, the cells -> tissues -> organs -> organ systems hierarchy, and
// the higher-level interactions between organ systems.
//
// Cell names, functions, tissues, organs, systems and system links are taken
// from the paper's lists. Germ-layer / lineage assignments use standard
// developmental biology; where the paper's grouping differs from the known
// developmental origin this is noted in `note` (e.g. osteoclasts and microglia
// are myeloid, Schwann cells and melanocytes come from the neural crest).

export const GERM_LAYERS = [
  { id: 'ectoderm', name: 'Ectoderm', color: 0x5aa8ff },
  { id: 'crest', name: 'Neural crest', color: 0x7bdff2 },
  { id: 'mesoderm', name: 'Mesoderm', color: 0xff6f91 },
  { id: 'endoderm', name: 'Endoderm', color: 0xffc36e },
  { id: 'germline', name: 'Germline', color: 0xc38bff },
];

const c = (name, fn = '', note = '') => ({ name, fn, note });

export const LINEAGES = [
  { id: 'blood', name: 'Hematopoietic (blood & immune)', layer: 'mesoderm', stem: 'Hematopoietic stem cell (HSC)', cells: [
    c('Erythrocyte (red blood cell)', 'Transport oxygen and carbon dioxide'), c('Neutrophil', 'Phagocytosis of pathogens'), c('Eosinophil', 'Granulocyte; anti-parasite defence'), c('Basophil', 'Granulocyte; histamine release'),
    c('Monocyte', 'Differentiates into macrophages and dendritic cells'), c('Macrophage', 'Phagocytosis, tissue repair'), c('Dendritic cell', 'Antigen presentation'),
    c('B cell', 'Adaptive immune response'), c('Plasma cell', 'Antibody-secreting B cell'), c('Memory B cell', 'Long-lasting immunity'),
    c('Helper T cell', 'Assist other white blood cells in immunologic processes'), c('Cytotoxic T cell', 'Kill infected cells'), c('Regulatory T cell', 'Suppress immune responses to maintain homeostasis'), c('Memory T cell', 'Rapid recall response'),
    c('Natural killer (NK) cell', 'Destroy virus-infected cells and tumours'), c('Megakaryocyte', 'Platelet-producing'), c('Platelet (thrombocyte)', 'Blood clotting'), c('Mast cell', 'Allergic and inflammatory responses'),
    c('Kupffer cell', 'Liver macrophages; phagocytose pathogens and debris'), c('Langerhans cell', 'Skin dendritic cells; present antigens to T cells'),
    c('Alveolar macrophage', 'Clear particles from alveoli'), c('Microglia', 'Resident immune cells of the CNS', 'Listed with glia in the paper; derived from yolk-sac myeloid progenitors'),
    c('Osteoclast', 'Bone resorption', 'Listed with bone cells in the paper; derived from the monocyte/macrophage lineage'),
  ] },
  { id: 'mesenchymal', name: 'Mesenchymal (connective tissue)', layer: 'mesoderm', stem: 'Mesenchymal stromal cell', cells: [
    c('White adipocyte', 'Fat storage'), c('Brown adipocyte', 'Thermogenesis'), c('Osteoblast', 'Bone formation'), c('Osteocyte', 'Mature bone cell; maintains bone'), c('Bone lining cell', 'Quiescent osteoblast on bone surfaces'),
    c('Periosteal cell', 'Bone growth and repair'), c('Chondrocyte', 'Cartilage formation'), c('Fibroblast', 'Secrete extracellular matrix and collagen'), c('Myofibroblast', 'Wound healing and fibrosis'),
    c('Tenocyte', 'Tendon cell'), c('Ligament cell', 'Maintain and repair ligament tissue'), c('Synovial cell', 'Produce synovial fluid for lubrication'), c('Bone marrow stromal cell', 'Support hematopoiesis'),
    c('Pericyte', 'Blood vessel support'), c('Hepatic stellate cell', 'Vitamin A storage, liver fibrosis'), c('Mesangial cell', 'Structural support in the glomerulus'), c('Follicular dendritic cell', 'Present antigen in lymph-node follicles', 'Stromal (not hematopoietic) origin'),
  ] },
  { id: 'muscle', name: 'Muscle', layer: 'mesoderm', stem: 'Satellite cell (muscle stem cell)', cells: [
    c('Skeletal muscle fibre', 'Contraction and movement'), c('Myoblast', 'Precursor that fuses into muscle fibres'), c('Satellite cell', 'Muscle repair and regeneration'),
    c('Vascular smooth muscle cell', 'Regulate vessel tone'), c('Intestinal smooth muscle cell', 'Peristalsis'), c('Uterine smooth muscle cell', 'Uterine contraction'), c('Bladder smooth muscle cell', 'Voiding'),
  ] },
  { id: 'cardiac', name: 'Cardiac', layer: 'mesoderm', stem: 'Cardiac progenitor', cells: [
    c('Cardiomyocyte', 'Heart muscle contraction'), c('Purkinje fibre', 'Specialised conduction'), c('Cardiac pacemaker cell', 'Set the heart rhythm'), c('Endocardial cell', 'Line the heart chambers'), c('Epicardial (mesothelial) cell', 'Outer heart layer'),
  ] },
  { id: 'endothelial', name: 'Endothelial & vascular', layer: 'mesoderm', stem: 'Hemangioblast / angioblast', cells: [
    c('Capillary endothelial cell', 'Exchange between blood and tissue'), c('Arterial endothelial cell', 'Line arteries'), c('Venous endothelial cell', 'Line veins'), c('Lymphatic endothelial cell', 'Line lymphatic vessels'),
    c('Liver sinusoidal endothelial cell', 'Fenestrated liver sinusoids'), c('Blood-brain barrier endothelial cell', 'Tight barrier of the CNS'), c('Glomerular endothelial cell', 'Fenestrated filtration barrier'),
  ] },
  { id: 'renal', name: 'Renal', layer: 'mesoderm', stem: 'Nephron progenitor (cap mesenchyme)', cells: [
    c('Podocyte', 'Kidney filtration'), c('Proximal tubule cell', 'Reabsorption'), c('Distal tubule cell', 'Ion balance'), c('Principal cell', 'Sodium and water balance'), c('Intercalated cell', 'Acid-base balance'),
    c('Juxtaglomerular cell', 'Secrete renin to regulate blood pressure'), c('Macula densa cell', 'Sense tubular NaCl'),
  ] },
  { id: 'gonadal', name: 'Reproductive (somatic)', layer: 'mesoderm', stem: 'Gonadal somatic progenitor', cells: [
    c('Leydig cell', 'Testosterone secretion'), c('Sertoli cell', 'Support and nourish developing sperm'), c('Granulosa cell', 'Support oocyte development'), c('Theca cell', 'Produce androgens'),
    c('Luteal cell', 'Produce progesterone after ovulation'), c('Cumulus cell', 'Surround the oocyte'), c('Decidual cell', 'Maternal part of the placenta'), c('Endometrial cell', 'Uterine lining'), c('Fallopian tube epithelial cell', 'Oocyte transport'),
  ] },
  { id: 'adrenalCortex', name: 'Adrenal cortex', layer: 'mesoderm', stem: 'Adrenocortical progenitor', cells: [
    c('Zona glomerulosa cell', 'Aldosterone'), c('Zona fasciculata cell', 'Cortisol'), c('Zona reticularis cell', 'Androgens'),
  ] },
  { id: 'neural', name: 'Neural (CNS)', layer: 'ectoderm', stem: 'Neural stem cell / radial glia', cells: [
    c('Sensory neuron', 'Transmit sensory information'), c('Motor neuron', 'Drive muscle contraction'), c('Interneuron', 'Connect neurons in circuits and reflexes'), c('Purkinje cell (cerebellum)', 'Motor control'),
    c('Granule cell (cerebellum)', 'Most numerous neuron; sensory processing'), c('Astrocyte', 'Support and protect neurons'), c('Oligodendrocyte', 'Myelinate CNS axons'), c('Ependymal cell', 'Line the ventricles and central canal'),
    c('Radial glial cell', 'Scaffold for neuron migration during development'), c('Magnocellular neurosecretory cell', 'Oxytocin / vasopressin'), c('Parvocellular neurosecretory cell', 'Hypothalamic releasing hormones'), c('Pinealocyte', 'Melatonin secretion'),
  ] },
  { id: 'crest', name: 'Neural crest', layer: 'crest', stem: 'Neural crest cell', cells: [
    c('Schwann cell', 'Myelinate peripheral axons'), c('Satellite glial cell', 'Surround neuron cell bodies in ganglia'), c('Enteric neuron', 'Gut nervous system'), c('Enteric glial cell', 'Support enteric neurons'),
    c('Dorsal root ganglion neuron', 'Somatosensory neuron'), c('Melanocyte', 'Pigment production'), c('Chromaffin cell', 'Epinephrine and norepinephrine'), c('Parafollicular (C) cell', 'Calcitonin'),
    c('Olfactory ensheathing cell', 'Guide olfactory axons'), c('Keratocyte (cornea)', 'Corneal stroma'),
  ] },
  { id: 'sensory', name: 'Sensory & neuroendocrine placodes', layer: 'ectoderm', stem: 'Placodal / retinal progenitor', cells: [
    c('Rod photoreceptor', 'Detect light'), c('Cone photoreceptor', 'Detect colour'), c('Retinal ganglion cell', 'Send visual signals to the brain'), c('Bipolar cell', 'Relay photoreceptor signals'), c('Amacrine cell', 'Retinal interneuron'), c('Horizontal cell', 'Lateral inhibition in the retina'),
    c('Inner hair cell (cochlea)', 'Detect sound'), c('Outer hair cell (cochlea)', 'Amplify sound'), c('Vestibular hair cell', 'Detect head position and motion'), c('Olfactory receptor neuron', 'Detect odour molecules'),
    c('Taste receptor cell', 'Sweet, sour, bitter, salty, umami'), c('Somatotroph', 'Growth hormone'), c('Lactotroph', 'Prolactin'), c('Corticotroph', 'ACTH'), c('Thyrotroph', 'TSH'), c('Gonadotroph', 'LH and FSH'), c('Melanotroph', 'MSH'),
  ] },
  { id: 'epidermal', name: 'Epidermal & glandular (surface ectoderm)', layer: 'ectoderm', stem: 'Epidermal (basal) stem cell', cells: [
    c('Keratinocyte', 'Form the outer layer of the skin'), c('Basal cell', 'Stem cells of the epithelium'), c('Corneocyte', 'Cornified surface cell'), c('Merkel cell', 'Detect light touch and texture'),
    c('Sebaceous gland cell', 'Secrete sebum'), c('Eccrine sweat gland cell', 'Water and small molecules'), c('Apocrine sweat gland cell', 'Odoriferous secretions'), c('Myoepithelial cell', 'Contract to expel glandular secretions'),
    c('Mammary gland cell', 'Secrete milk'), c('Lacrimal gland cell', 'Produce tears'), c('Ceruminous gland cell', 'Produce earwax'), c('Salivary serous cell', 'Secrete enzymes'), c('Salivary mucous cell', 'Secrete mucus'),
  ] },
  { id: 'gut', name: 'Gastrointestinal epithelium', layer: 'endoderm', stem: 'Intestinal (Lgr5+) / gastric stem cell', cells: [
    c('Enterocyte', 'Absorb nutrients'), c('Goblet cell', 'Secrete mucus'), c('Paneth cell', 'Secrete antimicrobial peptides'), c('Enteroendocrine G cell', 'Gastrin'), c('Enteroendocrine I cell', 'Cholecystokinin'),
    c('Enteroendocrine S cell', 'Secretin'), c('Enteroendocrine K cell', 'GIP'), c('Enterochromaffin-like cell', 'Histamine'), c('M cell', 'Antigen sampling over Peyer\'s patches'),
    c('Gastric chief cell', 'Pepsinogen'), c('Parietal cell', 'Hydrochloric acid'), c("Brunner's gland cell", 'Alkaline mucus in the duodenum'),
  ] },
  { id: 'hepatopancreatic', name: 'Liver & pancreas', layer: 'endoderm', stem: 'Hepatoblast / pancreatic progenitor', cells: [
    c('Hepatocyte', 'Detoxification, protein synthesis'), c('Cholangiocyte', 'Bile duct epithelium'), c('Acinar cell', 'Digestive enzymes'), c('Pancreatic ductal cell', 'Transport enzymes to the intestine'),
    c('Beta cell', 'Insulin'), c('Alpha cell', 'Glucagon'), c('Delta cell', 'Somatostatin'), c('PP cell', 'Pancreatic polypeptide'), c('Epsilon cell', 'Ghrelin'),
  ] },
  { id: 'respiratory', name: 'Respiratory epithelium', layer: 'endoderm', stem: 'Airway basal / club cell', cells: [
    c('Type I pneumocyte', 'Gas exchange'), c('Type II pneumocyte', 'Surfactant production'), c('Club cell', 'Protective secretions in bronchioles'), c('Ciliated epithelial cell', 'Move mucus out of the airways'),
    c('Airway goblet cell', 'Mucus to trap debris'), c('Tracheal gland cell', 'Airway secretions'),
  ] },
  { id: 'endocrine', name: 'Endocrine (endodermal)', layer: 'endoderm', stem: 'Pharyngeal endoderm progenitor', cells: [
    c('Thyroid follicular cell', 'Thyroxine'), c('Parathyroid chief cell', 'Parathyroid hormone (calcium)'), c('Thymic epithelial cell', 'T-cell education'), c('Prostate epithelial cell', 'Prostatic secretions'),
  ] },
  { id: 'germ', name: 'Germ cells', layer: 'germline', stem: 'Primordial germ cell / spermatogonial stem cell', cells: [
    c('Sperm cell', 'Male gamete'), c('Oocyte', 'Female gamete'), c('Spermatogonial stem cell', 'Self-renewing sperm precursor'),
  ] },
];

/** Morphology class used for the 3D glyph of each cell type. */
export function morphology(name) {
  const n = name.toLowerCase();
  if (/erythrocyte/.test(n)) return 'rbc';
  if (/platelet/.test(n)) return 'disc';
  if (/neutrophil|eosinophil|basophil/.test(n)) return 'lobed';
  if (/neuron|purkinje cell|granule cell|motor|interneuron|ganglion|bipolar|amacrine|horizontal|osteocyte|dendritic|astrocyte|oligodendro|microglia|melanocyte|podocyte|radial glial/.test(n)) return 'star';
  if (/muscle fibre|myoblast|cardiomyocyte|purkinje fibre|pacemaker/.test(n)) return 'capsule';
  if (/fibroblast|smooth muscle|tenocyte|ligament|pericyte|stellate|schwann|satellite|keratocyte|myoepithelial/.test(n)) return 'spindle';
  if (/adipocyte|megakaryocyte|oocyte/.test(n)) return 'big';
  if (/endothelial|type i pneumocyte|corneocyte|mesothelial|endocardial|squamous|bone lining/.test(n)) return 'flat';
  if (/enterocyte|goblet|columnar|ciliated|paneth|chief|parietal|hair cell|photoreceptor|club|tubule|principal|intercalated|keratinocyte|cholangiocyte|basal cell|acinar|ductal|follicular|pneumocyte/.test(n)) return 'columnar';
  if (/sperm/.test(n)) return 'sperm';
  return 'round';
}

export const TISSUES = [
  { name: 'Epithelial tissue', kinds: ['Simple squamous', 'Simple cuboidal', 'Simple columnar', 'Stratified squamous', 'Transitional'], cells: ['Keratinocyte', 'Enterocyte', 'Goblet cell', 'Paneth cell', 'Type I pneumocyte', 'Proximal tubule cell'] },
  { name: 'Connective tissue', kinds: ['Loose (areolar)', 'Dense', 'Adipose', 'Reticular', 'Cartilage', 'Bone', 'Blood'], cells: ['Fibroblast', 'White adipocyte', 'Chondrocyte', 'Osteocyte', 'Erythrocyte (red blood cell)', 'Macrophage', 'Mast cell'] },
  { name: 'Muscle tissue', kinds: ['Skeletal', 'Cardiac', 'Smooth'], cells: ['Skeletal muscle fibre', 'Cardiomyocyte', 'Vascular smooth muscle cell'] },
  { name: 'Nervous tissue', kinds: ['Neurons', 'Neuroglia'], cells: ['Sensory neuron', 'Motor neuron', 'Astrocyte', 'Oligodendrocyte', 'Microglia', 'Ependymal cell'] },
];

export const ORGANS = [
  { name: 'Skin', tissues: ['Epithelial tissue', 'Connective tissue'], parts: 'Epidermis (keratinocytes, melanocytes, Langerhans cells); dermis (fibroblasts, macrophages, mast cells); hypodermis (adipocytes)', system: 'Integumentary System' },
  { name: 'Heart', tissues: ['Muscle tissue', 'Epithelial tissue'], parts: 'Myocardium (cardiac muscle cells), endocardium (endothelial cells), epicardium (mesothelial cells)', system: 'Cardiovascular System' },
  { name: 'Liver', tissues: ['Epithelial tissue', 'Connective tissue'], parts: 'Hepatocytes, Kupffer cells, stellate cells, bile duct cells', system: 'Digestive System' },
  { name: 'Kidney', tissues: ['Epithelial tissue', 'Connective tissue'], parts: 'Nephrons (podocytes, proximal/distal tubule, collecting duct cells); renal corpuscle (glomerular endothelial, mesangial cells)', system: 'Urinary System' },
  { name: 'Lungs', tissues: ['Epithelial tissue', 'Connective tissue'], parts: 'Alveoli (type I/II pneumocytes, alveolar macrophages); bronchioles (ciliated, goblet cells)', system: 'Respiratory System' },
  { name: 'Pancreas', tissues: ['Epithelial tissue'], parts: 'Islets of Langerhans (alpha, beta, delta, PP cells); acinar cells; ductal cells', system: 'Endocrine System' },
  { name: 'Stomach', tissues: ['Epithelial tissue', 'Connective tissue', 'Muscle tissue'], parts: 'Mucosa (chief, parietal, goblet, enteroendocrine cells); submucosa; muscularis (smooth muscle)', system: 'Digestive System' },
  { name: 'Intestines', tissues: ['Epithelial tissue', 'Connective tissue', 'Muscle tissue'], parts: 'Mucosa (enterocytes, goblet, Paneth, enteroendocrine cells); submucosa; muscularis', system: 'Digestive System' },
  { name: 'Spleen', tissues: ['Connective tissue'], parts: 'White pulp (lymphocytes); red pulp (macrophages, red blood cells)', system: 'Lymphatic/Immune System' },
  { name: 'Brain', tissues: ['Nervous tissue'], parts: 'Neurons and neuroglia', system: 'Nervous System' },
];

export const SYSTEMS = [
  { name: 'Integumentary System', organs: 'Skin, hair, nails, glands (sebaceous, sweat)' },
  { name: 'Muscular System', organs: 'Skeletal muscles, cardiac muscle, smooth muscle' },
  { name: 'Nervous System', organs: 'Brain, spinal cord, peripheral nerves' },
  { name: 'Endocrine System', organs: 'Pituitary, thyroid, adrenal glands, pancreas' },
  { name: 'Cardiovascular System', organs: 'Heart, blood vessels (arteries, veins, capillaries)' },
  { name: 'Respiratory System', organs: 'Lungs, trachea, bronchi, alveoli' },
  { name: 'Digestive System', organs: 'Mouth, esophagus, stomach, intestines, liver, pancreas' },
  { name: 'Urinary System', organs: 'Kidneys, ureters, bladder, urethra' },
  { name: 'Reproductive System', organs: 'Testes, ovaries, uterus, prostate, mammary glands' },
  { name: 'Lymphatic/Immune System', organs: 'Lymph nodes, spleen, thymus, bone marrow' },
];

// the paper's higher_level_connections with its descriptions
export const SYSTEM_LINKS = [
  { a: 'Cardiovascular System', b: 'Respiratory System', text: 'The heart pumps oxygenated blood received from the lungs to the body; deoxygenated blood returns to the lungs.' },
  { a: 'Nervous System', b: 'Muscular System', text: 'Neurons transmit signals to muscle fibres to initiate movement and coordination.' },
  { a: 'Endocrine System', b: 'Digestive System', text: 'The pancreas secretes insulin and glucagon; hormones regulate digestive processes.' },
  { a: 'Lymphatic/Immune System', b: 'Cardiovascular System', text: 'Immune cells travel through lymphatic vessels and blood to sites of infection; lymph nodes filter lymph.' },
  { a: 'Urinary System', b: 'Cardiovascular System', text: 'The kidneys filter blood to remove waste and maintain fluid balance; blood supplies the kidneys for filtration.' },
  { a: 'Reproductive System', b: 'Endocrine System', text: 'Hormones regulate reproductive processes; gonads produce sex hormones that act on other systems.' },
  { a: 'Nervous System', b: 'Endocrine System', text: 'The hypothalamus regulates hormone release from the pituitary gland.' },
];

export function allCellTypes() {
  return LINEAGES.flatMap((l) => l.cells.map((cell) => ({ ...cell, lineage: l.id, layer: l.layer, stem: l.stem })));
}
