# UTC AI-Powered Whole-Body Simulator: Upgrade Blueprint

**Status:** proposed architecture, not yet implemented. **Repository:** Mecca-Research/StemCellSimulator. **Baseline:** main, November 2026.

## Architecture decision

Preserve the deterministic Three.js/JavaScript simulator and extend it with TWO tracks:

1. **Biological Knowledge Infrastructure (BKI)**: federated scientific data access, source-pinned snapshots, biological entity normalization, provenance-aware knowledge graph, and reliable registry releases.
2. **Hybrid Simulation Infrastructure (HSI)**: versioned deterministic model runs, external scientific simulators, a server-side AI model broker, uncertainty verification and data-driven 3D animations.

The existing GitHub Pages demo must continue working offline. AI-generated mechanisms remain hypotheses until supported by evidence; the architecture must not claim that simulations are validated substitutes for clinical trials.

## Current repository and required repairs

Already usable:
- js/models/core/ode.js and rng.js: ODE, stochastic numerical methods and seeded randomness.
- js/models/pathways, morpho, protein and organ: pure simulation modules separated from rendering.
- js/models/network/registry.js: module-defined reaction lists, bipartite graph and PCA display.
- js/models/atlas/atlas.js: existing presentation ontology, with 164 named cell types.
- js/scenes and js/engine: Three.js display, microscopy-style modes, controls and scene lifecycle.
- tools/fetch_micrographs.py and tools/build_assets.py: SHA-256-pinned evidence downloads and reproducible visualization assets.
- tests/models, tests/tools and tests/e2e: foundational regression tests.

Required:
- The graph currently keys species by trimmed **display name**, causing false merges across tissues, chemical states, organisms and subcellular compartments. Replace with stable, qualified entity IDs and context-specific species instances.
- Module-declared REACTIONS, experimentally curated biochemical reactions, and 2024 paper estimates are DIFFERENT kinds of evidence. Keep them separate in the registry and visualization.
- The 2024 paper's reported 500/1,051/2,051/3,021 reactions arose in part from synthetic example code; mark them as illustrative, not experimentally validated counts.
- The hard-coded atlas is not an authenticated Human Cell Atlas release.
- Validate stiff kinetic regimes with an appropriate solver; Dormand-Prince RK5(4) is adaptive but not generally stiff-stable.

## 1. Database acquisition: choose the access mode by dataset

| Category | Initial providers | Connection mode |
|---|---|---|
| Curated biochemical reactions | Rhea, Reactome, ChEBI | Versioned bulk exports; selective REST/SPARQL lookups |
| Protein and gene IDs | UniProt, Ensembl, HGNC | Pinned bulk tables plus targeted API lookup |
| Cell identity and ontology | Cell Ontology, UBERON, GO via OLS4 | Versioned ontology snapshot and API |
| Human single-cell expression | HCA, CELLxGENE Census | Cloud-native TileDB-SOMA slices and pinned subset exports, not full local mirror |
| Spatial and organ context | HuBMAP, Human Protein Atlas, Allen | Metadata APIs and selective licensed files |
| Chemical and bioactivity data | PubChem, ChEMBL, GPCRdb | Bulk assay snapshots plus cached REST lookups |
| Dynamic systems models | BioModels/SBML, Physiome/CellML | Download versioned models for scientific runners |
| Perturbations and safety | LINCS, JUMP, EPA CompTox | Snapshot releases and assay-filtered extracts |

A proper connector speaks each provider's supported interface: REST, SPARQL, Python client, bulk export, signed object-store downloads, or standards-based files. API calls are not the only option.

Never put indiscriminate scientific data downloads in the browser. Public APIs are shared services and may have quotas, version-specific schemas, and different licenses. PubChem asks clients to stay within a published 5-requests-per-second ceiling; favor bulk downloads for high-volume work.

Canonical ingestion workflow:
Discover -> Fetch -> Verify Checksum/License -> Parse -> Normalize -> Crosswalk -> Quality Gate -> Review Conflicts -> Publish Immutable Dataset Release.

Use a manifest recording provider, endpoint, source release, timestamp, content hash, license, source attribution and processing version. No provider outage may silently trigger synthetic fallback data.

## 2. Storage architecture: cloud and local, not GitHub as a database

- **GitHub:** source code, schemas, tests, small licensed fixtures, manifests and model cards. No giant single-cell matrices, API secrets or private donor data.
- **S3-compatible object storage:** original immutable source releases; extracted Parquet/AnnData; models; run reports; simulation trajectories; hashed asset bundles.
- **PostgreSQL (canonical registry):** stable entity IDs, synonyms, ontology mappings, source records, evidence assertions, context, parameter sets, datasets, model versions, simulation run metadata.
- **DuckDB + Parquet:** analytical exploration, distributed-ready columnar processing and source-release comparisons.
- **Optional graph database:** only introduce Neo4j or similar after relational typed edges and graph projections prove insufficient for workload.
- **Optional vector search:** paper/model retrieval, never a replacement for authoritative source references.

Begin with a small FastAPI service + Postgres + object storage and a background worker. Avoid early infrastructure sprawl.

## 3. Identity and context reconciliation

Define typed entities: compound, compound-state, gene, protein, isoform, protein-state, biological reaction, signaling relation, cell type, cell instance, tissue, organ, experiment, assay, observation, dataset, model, run, source and hypothesis.

Use domain identifiers:
- Compounds: ChEBI/PubChem/ChEMBL IDs plus InChI where appropriate; preserve stereo/charge/tautomer/isotope/salt differences.
- Genes: Ensembl/HGNC and explicit NCBI taxonomy.
- Proteins: UniProt with isoform and post-translational-state representation.
- Cell identities: Cell Ontology with species, tissue, development and experimental state; anatomy from UBERON, compartments from GO.
- Reactions: Rhea/Reactome provenance, explicit stoichiometry, direction, catalysts, modifiers, location, species and rates when known.
- Literature and experiments: DOI/PMID, assay, method, organism, measured effect, error, detection limit, units and source version.

**Never merge records just because labels look alike.** Match with identifier crosswalks and evidence; preserve unresolved equivalence as possible_match. Keep a chemical's molecular identity separate from the same chemical's **distinct cytosolic, membrane and extracellular state variables**.

Chemical reactions are hyperedges with all participants and stoichiometries. A receptor association, gene regulation and ligand binding observation are different typed relations, not all equivalent "reactions."

## 4. Evidence registry and confidence

Every assertion stores: claim ID; subject/predicate/object; numerical values with units; source accession and dataset version; assay context; species, cell/tissue/compartment; measured condition; time; evidence category; sample count and error; reviewer; license; model/version if predicted.

Evidence classes: measured, experimentally curated, literature supported, mechanistic inference, statistical prediction, illustrative, contradicted, rejected.

Keep conflicting results as separate claims. Do not collapse contradiction into a fabricated consensus. Support claim-scoped approval states; only reviewed evidence can promote a hypothesis into the trusted graph.

Use separate metrics for **evidence strength** and **predictive reliability**. Probability values must come from validated calibrated statistical models; an LLM confidence number is not a scientific confidence interval.

Quality checks: semantic ID collisions, source checksum, reference integrity, duplicate assay records, unit normalization, mass/charge and stoichiometric balance where applicable, species/compartment consistency, ontology validity, licensing and release drift.

## 5. Hybrid deterministic simulation

Keep existing plain-JS simulation modules. Add a versioned, testable transport contract:

**RunSpec:** run ID, knowledge release, model ID/version/hash, biological context, initial state, parameter values+units, interventions, solver/version, random seed, simulation time and requested outputs.

**RunResult:** trajectories with units and timestamps, state checkpoints, model warnings, numerical error/conservation tests, optional uncertainty intervals, hashes, source/parameter provenance and scientific validity status.

**FrameState:** time, object identities, spatial position/morphology, visible signals and links back to observed or predicted values. Renderer consumes FrameState, not LLM prose.

Run deterministic local simulation first. Introduce an API worker for large runs and standardized model interoperability (SBML/COMBINE, CellML) later. Evaluate COPASI, VCell, BioModels and PK-Sim/MoBi as external engines when the relevant capability is validated. Explicitly separate mathematical correctness from biological fidelity.

## 6. AI model hosting and orchestration

Architecture:
Browser -> authenticated backend -> evidence retrieval -> model broker -> typed, cited proposal -> verifier -> report / reviewed simulation request.

- Public GitHub Pages must NEVER contain AI provider keys.
- Provider-independent model adapter with model/version, rate limits, timeouts, token/cost quotas, request hashing, response caching and audit logs.
- AI tools limited initially to source-backed read queries, graph walks, scenario explanation, proposing hypotheses, and invoking approved simulation jobs.
- Require all claims to reference exact registry IDs. Missing evidence creates a flagged gap, not an invented edge or parameter.
- LLM predictions cannot directly rewrite approved kinetics, model parameters, canonical records or existing experimental observations.
- Specialized statistical/chemical/cellular models compute numerical predictions; LLMs coordinate, retrieve and explain. Server-side review gate controls promotion.
- Protect private human data with consent and lawful processing, authorization, encryption and controlled external provider sharing.

## 7. Scientific visualization upgrade

Preserve 11 existing scene entries, 3D/confocal/H&E modes and offline behavior. Add:
- Evidence inspector with source IDs, organism, tissue, measured vs inferred vs illustrative status and uncertainty.
- Queryable reaction explorer with typed, filterable subgraphs and evidence conflict overlays.
- Cell and organ atlas with linked ontology IDs and spatial constraints.
- Run timeline showing modeled state trajectories, concentration/kinetic charts, parameter changes and 3D animations in sync.
- Observed vs predicted overlays, simulation ensembles and uncertainty bands.
- Playback provenance: a static microscopy image is evidence for geometry, not evidence of its future motion. A reconstructed animation must be labeled modeled/illustrative unless supported by time-lapse data.
- Server-side graph filtering and progressive rendering for large graphs; do not reuse current dense PCA approach on millions of nodes.

## 8. Model validation and ethical requirements

Validation ladder: (1) unit and conservation tests; (2) curated known-model reproduction; (3) data-source- and lab-separated retrospective benchmarks; (4) independent prospective human-cell predictions; (5) tissue and organ specific assessment; (6) regulated human studies where justified.

Publish per-endpoint MAE/RMSE, classification/calibration metrics, interval coverage, out-of-distribution warnings, error analysis, model version and reproducibility results. No universal "99% accuracy" claim.

Research workflow prioritizes well-validated nonanimal experiments; any proposed animal trial requires explicit scientific necessity, alternatives assessment, independent ethics review and welfare refinement. AI simulations are not alone sufficient to approve a clinical trial.

## 9. First vertical slice and acceptance criteria

Use existing Wnt/Notch and reaction network, then connect a *tiny version-pinned* Rhea + Reactome + UniProt + Cell Ontology fixture.

Must demonstrate:
1. One live scientific source can be imported reproducibly with license/checksum/version.
2. Two unrelated source records can be matched through stable IDs without accidentally merging compartmentalized species.
3. An unresolved cross-source contradiction survives as two auditable claims.
4. Curated measured reactions are shown visually distinct from code-defined phenomenological reactions and 2024 illustrative paper counts.
5. A deterministic local RunSpec can be reproduced; recorded outputs and tolerances pass CI.
6. A server-side AI response explains a simulated change with source IDs; fabricated claims are quarantined.
7. A run report contains provenance, parameter/solver versions, known limitations and uncertainty.
8. All existing scenes continue to render, including offline GitHub Pages.

Only then extend the verified graph to known indole/serotonin pharmacology, and later integrate tissue exposure, neuronal and immune dynamics.

## 10. Suggested additive repo layout

- schemas/: JSON Schemas and RunSpec/RunResult/FrameState, unit and ontology conventions.
- backend/app/: FastAPI query, evidence, simulation, run and AI endpoints.
- backend/ingestion/: source connectors, version pinning, ETL and crosswalks.
- backend/registry/: SQL migrations, evidence lifecycle, conflict resolution.
- backend/simulation/: deterministic and third-party model adapters.
- backend/ai/: provider broker, tool access, typed result verification.
- backend/verification/: benchmarks, calibration, source auditing.
- js/client/: optional remote API client, playback adapter, provenance panels.
- data/manifests/ and data/fixtures/: small test datasets, not giant corpora.
- tests/: scientific contracts, data tests, API security, CI regressions, existing E2E.

## 11. Phased PR plan

**PR A, foundation:** typed entities/reactions, RunSpec schemas, provenance labels, graph identity migration. No backend required.

**PR B, source connectors:** Rhea/Reactome imports, source manifests, Postgres migrations, conflict-aware crosswalks, fixture tests.

**PR C, cell atlas:** ontology and UniProt linking, version-pinned selective HCA/CELLxGENE Census imports, no full-dataset mirroring.

**PR D, simulation service:** deterministic job API, cached results, storage, replay and numerical tests.

**PR E, AI broker:** hosted provider adapter, structured citations, hypothesis quarantine, validation and security tests.

**PR F, UI:** evidence inspector, source-backed graph, synced timeline, multiple scenario playback.

**PR G, indole case study:** measured receptor selectivity and signaling benchmarks before any claim of novel molecule efficacy.

**Later:** explicitly validated organ coupling, PBPK/external model adapters, and prospective independent research collaboration.

## 12. Official upstream programmatic references

- Rhea: https://www.rhea-db.org/help/api
- Reactome: https://reactome.org/download/
- CELLxGENE Census: https://chanzuckerberg.github.io/cellxgene-census/python-api.html
- Ontology Lookup Service 4: https://www.ebi.ac.uk/ols4/about
- ChEMBL: https://www.ebi.ac.uk/chembl/api/data/docs
- PubChem: https://pubchem.ncbi.nlm.nih.gov/docs/programmatic-access
- Open Systems Pharmacology: https://www.open-systems-pharmacology.org/
- BioModels: https://www.ebi.ac.uk/biomodels/

This is a reviewed architecture **plan**, not a claim that new database connectors, cloud runners, or clinical validation have been implemented.
